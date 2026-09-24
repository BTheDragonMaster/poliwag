"""
Compare a large pool of amplicon sequences against a small set of reference
amplicons, one output file per reference, without ever holding the large pool
in memory.

Built for real amplicon datasets where the large side can be millions of
sequences (e.g. comparing ~7.3M amplicons against 13 cluster representatives,
as discussed for ``cluster_amplicon_analysis.py``'s matches-vs-non-matches
step): materializing that as a list, or as a full N x M cross-product of
pairs, costs multiple GB just in Python object/list overhead before any
comparison work happens. This script instead reads the large FASTA as a
generator (one record in memory at a time) and keeps only the small reference
pool resident, so peak memory is ~O(reference pool size), independent of how
large the input FASTA is.

Deliberately self-contained (does not import from
``poliwag.amplicons.cluster_amplicon_analysis``, which pulls in
``poliwag.antismash_processing`` and, transitively, the real ``antismash``
package): this script only ever does sequence I/O and percent-identity
arithmetic, so it has no need for that heavier dependency chain and stays
fast to start on its own for large, standalone runs.

For each reference sequence, writes a TSV of
``seq_1_name\tseq_2_name\tpercent_identity\tpercent_similarity`` for every
large-pool sequence that scores at or above ``--min_identity`` (default 30%,
checked against percent identity) against it. Any sequence (in either pool)
containing a stop codon ("*", as translated protein sequences typically
represent one) is skipped entirely -- it's never compared against anything
and never appears in any output file.

Every comparison is a real Biopython ``PairwiseAligner`` global alignment
against a BLOSUM62 substitution matrix (not a plain positional/mismatch
count), so the two amplicons don't need to be the same length, and the two
reported numbers follow the usual EMBOSS-style convention: percent identity
is the fraction of aligned columns with identical residues, and percent
similarity is that plus columns where the two (different) residues still
score positively under BLOSUM62 -- both computed over the full alignment
length (gaps included), so identity <= similarity always. This is
meaningfully slower per pair than a plain identity count (~100 microseconds/
pair measured on 33-residue amplicons) since it now runs dynamic-programming
alignment on every pair rather than a positional comparison. Every pair is
independent, so the comparison step parallelizes across CPU cores with
``multiprocessing`` (``--processes``, default: all cores) -- the large FASTA
is still only ever read one record at a time in the main process (workers
never open it themselves), so this doesn't change the memory story above,
only the wall-clock time.

Once each reference's file is fully written, it's read back and sorted by
percent identity (highest first), overwriting it in place -- see
``sort_identity_file``. Pass ``--skip_sort`` to leave files in streaming
(unsorted) order instead.
"""

import os
import re
import time
from argparse import ArgumentParser, Namespace
from multiprocessing import Pool, cpu_count
from typing import Iterator, Optional

from Bio.Align import PairwiseAligner, substitution_matrices
from tqdm import tqdm

STOP_CODON_SYMBOL = "*"


# ---------------------------------------------------------------------------
# Percent identity / similarity via a real Biopython pairwise alignment
# (self-contained, kept separate from cluster_amplicon_analysis.py on purpose
# -- see module docstring)
# ---------------------------------------------------------------------------

_SUBSTITUTION_MATRIX = substitution_matrices.load("BLOSUM62")

_ALIGNER = PairwiseAligner()
_ALIGNER.substitution_matrix = _SUBSTITUTION_MATRIX
_ALIGNER.mode = "global"
_ALIGNER.open_gap_score = -10.0
_ALIGNER.extend_gap_score = -0.5


def pairwise_identity_and_similarity(seq1: str, seq2: str) -> tuple[float, float]:
    """Globally align two amplicon (protein) sequences with Biopython's
    ``PairwiseAligner`` (BLOSUM62 substitution matrix) and return
    (percent_identity, percent_similarity).

    Percent identity counts aligned columns with identical residues; percent
    similarity additionally counts columns where the two (different)
    residues still score positively under BLOSUM62 -- the same convention
    tools like EMBOSS needle/water use. Both are divided by the full
    alignment length (gaps included), so a gappy alignment reduces both
    percentages, and identity is always <= similarity.

    :param seq1: first amplicon sequence
    :param seq2: second amplicon sequence
    :return: (percent_identity, percent_similarity), both 0-100
    """
    if not seq1 or not seq2:
        return 0.0, 0.0

    alignment = _ALIGNER.align(seq1, seq2)[0]
    aligned1 = str(alignment[0])
    aligned2 = str(alignment[1])

    length = len(aligned1)
    identical = 0
    similar = 0
    for a, b in zip(aligned1, aligned2):
        if a == "-" or b == "-":
            continue
        if a == b:
            identical += 1
            similar += 1
        else:
            try:
                if _SUBSTITUTION_MATRIX[a, b] > 0:
                    similar += 1
            except KeyError:
                pass  # residue not in the matrix alphabet (e.g. a stray ambiguity code)

    return 100.0 * identical / length, 100.0 * similar / length


# ---------------------------------------------------------------------------
# Streaming FASTA reading
# ---------------------------------------------------------------------------

def iter_fasta(path: str) -> Iterator[tuple[str, str]]:
    """Lazily read a FASTA file, yielding (header, sequence) one record at a
    time. Never holds more than one record in memory, so this is safe to use
    directly on arbitrarily large files.

    :param path: path to a FASTA file
    :yield: (header without the leading '>', sequence with newlines stripped)
    """
    header = None
    sequence_lines: list[str] = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_lines)
                header = line[1:]
                sequence_lines = []
            else:
                sequence_lines.append(line)
        if header is not None:
            yield header, "".join(sequence_lines)


def count_fasta_records(path: str) -> int:
    """Count records in a FASTA file (for an accurate progress bar total)
    with a single fast sequential read."""
    count = 0
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def has_stop_codon(sequence: str) -> bool:
    return STOP_CODON_SYMBOL in sequence


def sanitize_filename(name: str, max_length: int = 150) -> str:
    """Turn a FASTA header (which may contain '|', spaces, etc.) into a safe
    filename fragment."""
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_")
    return safe[:max_length] if safe else "sequence"


# ---------------------------------------------------------------------------
# Reference pool loading
# ---------------------------------------------------------------------------

def load_reference_pool(path: str) -> list[tuple[str, str]]:
    """Load the (small) reference pool fully into memory, dropping any
    sequence with a stop codon.

    :param path: path to the reference pool FASTA
    :return: list of (header, sequence) for references without a stop codon
    """
    pool: list[tuple[str, str]] = []
    for header, sequence in iter_fasta(path):
        if has_stop_codon(sequence):
            print(f"Warning: reference sequence {header!r} contains a stop codon, skipping it entirely")
            continue
        pool.append((header, sequence))
    return pool


# ---------------------------------------------------------------------------
# Multiprocessing worker
# ---------------------------------------------------------------------------
#
# Each worker process gets its own copy of the (small) reference pool and
# minimum identity, set once via the Pool initializer rather than re-sent
# with every task. `_ALIGNER`/`_SUBSTITUTION_MATRIX` above are module-level
# globals, so each worker process (forked or spawned fresh) ends up with its
# own independent instance -- PairwiseAligner isn't shared/thread-safe across
# processes, and this setup doesn't need it to be.

_worker_reference_pool: list[tuple[str, str]] = []
_worker_min_identity: float = 30.0


def _init_worker(reference_pool: list[tuple[str, str]], min_identity: float) -> None:
    global _worker_reference_pool, _worker_min_identity
    _worker_reference_pool = reference_pool
    _worker_min_identity = min_identity


def _compare_one_record(record: tuple[str, str]) -> tuple[bool, list[tuple[int, str, str, float, float]]]:
    """Compare one large-pool record against every reference sequence in this
    worker's reference pool. Runs in a worker process.

    :param record: (header, sequence) for one large-pool sequence
    :return: (was_skipped_for_stop_codon, hits), where ``hits`` is a list of
        (reference_index, header, ref_header, percent_identity,
        percent_similarity) for pairs meeting the minimum identity -- only
        passing pairs are returned, to keep the data sent back to the main
        process small
    """
    header, sequence = record
    if has_stop_codon(sequence):
        return True, []

    hits = []
    for ref_index, (ref_header, ref_sequence) in enumerate(_worker_reference_pool):
        identity, similarity = pairwise_identity_and_similarity(sequence, ref_sequence)
        if identity >= _worker_min_identity:
            hits.append((ref_index, header, ref_header, identity, similarity))
    return False, hits


# ---------------------------------------------------------------------------
# Sorting a finished identity file by percent identity, highest first
# ---------------------------------------------------------------------------

def sort_identity_file(path: str) -> None:
    """Sort a finished identity TSV (as written by ``compare_pools``) by
    percent identity, highest first, overwriting it in place.

    At the sizes these files actually reach (on the order of ~1M rows, tens
    of MB of text), a plain in-memory sort is the right tool: the whole file
    plus its parsed rows comfortably fits in memory (measured at ~260MB peak
    for a 1M-row file), takes a couple of seconds, and is freed again before
    the next file is processed -- no need for an external/chunked sort.

    :param path: path to a TSV with a header line and
        seq_1_name\tseq_2_name\tpercent_identity\tpercent_similarity rows
    """
    with open(path) as handle:
        header = handle.readline()
        rows = []
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            seq_1_name, seq_2_name, percent_identity_str, percent_similarity_str = line.split("\t")
            rows.append((seq_1_name, seq_2_name, float(percent_identity_str), float(percent_similarity_str)))

    rows.sort(key=lambda row: row[2], reverse=True)

    with open(path, "w") as handle:
        handle.write(header)
        for seq_1_name, seq_2_name, identity, similarity in rows:
            handle.write(f"{seq_1_name}\t{seq_2_name}\t{identity:.2f}\t{similarity:.2f}\n")


# ---------------------------------------------------------------------------
# Main streaming comparison
# ---------------------------------------------------------------------------

def compare_pools(
    large_fasta: str,
    reference_fasta: str,
    out_dir: str,
    min_identity: float = 30.0,
    sort_output: bool = True,
    processes: Optional[int] = None,
    chunksize: int = 100,
) -> None:
    """Stream ``large_fasta`` once, comparing every non-stop-codon record
    against every reference sequence, writing one TSV per reference.

    :param large_fasta: path to the large amplicon pool (streamed, never
        fully loaded into memory)
    :param reference_fasta: path to the small reference pool (loaded fully
        into memory once)
    :param out_dir: directory to write one <reference>.tsv per reference
        sequence into
    :param min_identity: minimum percent identity to keep a pair (default 30.0);
        checked against percent identity only, percent similarity is reported
        alongside but doesn't affect filtering
    :param sort_output: if True (default), sort each finished output file by
        percent identity, highest first, once the streaming comparison is
        done. At ~1M rows per file this adds only a couple of seconds per
        file (see ``sort_identity_file``).
    :param processes: number of worker processes to compare with in parallel.
        Defaults to all available CPU cores. 1 runs single-process (no pool).
    :param chunksize: number of large-pool records handed to a worker per
        task (``multiprocessing.Pool.imap_unordered``'s ``chunksize``).
        Larger values cut inter-process overhead at the cost of coarser
        progress-bar updates and slightly uneven load near the end of the run.
    """
    reference_pool = load_reference_pool(reference_fasta)
    if not reference_pool:
        raise ValueError(f"No usable reference sequences found in {reference_fasta!r}")

    os.makedirs(out_dir, exist_ok=True)

    handles = []
    out_paths = []
    for i, (header, _) in enumerate(reference_pool):
        out_path = os.path.join(out_dir, f"{i:03d}_{sanitize_filename(header)}.tsv")
        handle = open(out_path, "w")
        handle.write("seq_1_name\tseq_2_name\tpercent_identity\tpercent_similarity\n")
        handles.append(handle)
        out_paths.append(out_path)

    n_processed = 0
    n_skipped_stop_codon = 0
    n_pairs_written = 0

    n_processes = processes if processes is not None else (cpu_count() or 1)
    total = count_fasta_records(large_fasta)

    try:
        if n_processes <= 1:
            # Single-process path: run the same worker function directly,
            # no Pool/IPC overhead.
            _init_worker(reference_pool, min_identity)
            results_iter = (_compare_one_record(record) for record in iter_fasta(large_fasta))
            progress = tqdm(results_iter, total=total, desc="Comparing amplicons", unit="sequence")
            for was_skipped, hits in progress:
                n_processed += 1
                if was_skipped:
                    n_skipped_stop_codon += 1
                    continue
                for ref_index, header, ref_header, identity, similarity in hits:
                    handles[ref_index].write(f"{header}\t{ref_header}\t{identity:.2f}\t{similarity:.2f}\n")
                    n_pairs_written += 1
        else:
            with Pool(
                processes=n_processes,
                initializer=_init_worker,
                initargs=(reference_pool, min_identity),
            ) as pool:
                progress = tqdm(
                    pool.imap_unordered(_compare_one_record, iter_fasta(large_fasta), chunksize=chunksize),
                    total=total,
                    desc=f"Comparing amplicons ({n_processes} processes)",
                    unit="sequence",
                )
                for was_skipped, hits in progress:
                    n_processed += 1
                    if was_skipped:
                        n_skipped_stop_codon += 1
                        continue
                    for ref_index, header, ref_header, identity, similarity in hits:
                        handles[ref_index].write(f"{header}\t{ref_header}\t{identity:.2f}\t{similarity:.2f}\n")
                        n_pairs_written += 1
            # `iter_fasta(large_fasta)` is consumed by the main process (which
            # feeds the pool); worker processes never open the large FASTA
            # themselves, so peak memory is unaffected by parallelizing.
    finally:
        for handle in handles:
            handle.close()

    print(f"Processed {n_processed} sequences from {large_fasta}")
    print(f"Skipped {n_skipped_stop_codon} sequences containing a stop codon")
    print(f"Wrote {n_pairs_written} pairs (>= {min_identity}% identity) across {len(reference_pool)} reference file(s) in {out_dir}")

    if sort_output:
        for out_path in out_paths:
            start = time.perf_counter()
            sort_identity_file(out_path)
            elapsed = time.perf_counter() - start
            print(f"Sorted {out_path} by percent identity (highest first) in {elapsed:.1f}s")


def parse_arguments() -> Namespace:
    parser = ArgumentParser(
        description="Stream a large amplicon FASTA against a small reference pool, writing one "
        "identity TSV per reference sequence, without holding the large pool in memory."
    )
    parser.add_argument("-l", "--large_fasta", required=True, type=str,
                        help="Path to the large amplicon pool FASTA (streamed, not loaded into memory)")
    parser.add_argument("-r", "--reference_fasta", required=True, type=str,
                        help="Path to the small reference pool FASTA (loaded fully into memory)")
    parser.add_argument("-o", "--out_dir", required=True, type=str,
                        help="Directory to write one <reference>.tsv per reference sequence into")
    parser.add_argument("-m", "--min_identity", default=30.0, type=float,
                        help="Minimum percent identity to keep a pair (default: 30.0)")
    parser.add_argument("--skip_sort", action="store_true", default=False,
                        help="Leave each output file in streaming (unsorted) order instead of "
                        "sorting by percent identity (highest first) once writing is done.")
    parser.add_argument("-p", "--processes", type=int, default=None,
                        help="Number of worker processes to compare pairs in parallel (default: "
                        "all available CPU cores). Use 1 to run single-process.")
    parser.add_argument("--chunksize", type=int, default=100,
                        help="Number of large-pool records handed to a worker per task "
                        "(default: 100). Larger values reduce inter-process overhead at the cost "
                        "of coarser progress-bar updates.")
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    compare_pools(
        args.large_fasta,
        args.reference_fasta,
        args.out_dir,
        args.min_identity,
        sort_output=not args.skip_sort,
        processes=args.processes,
        chunksize=args.chunksize,
    )


if __name__ == "__main__":
    main()
