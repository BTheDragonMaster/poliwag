"""
Locate the region in each protein sequence that best matches a degenerate
DNA primer, by locally aligning a single consensus peptide translated from
the primer against every input sequence, and extract a fixed-length
amplicon anchored to where the *full* consensus peptide would start in the
target.

``degenerate_primer_amplicons.py`` requires an exact match to one of the
amino acid combinations the primer's degenerate codons encode. This script
instead tolerates mismatches: for each of the primer's codons it picks the
single most common amino acid the codon encodes (ties broken
alphabetically) to build one representative "consensus" peptide, then
finds the best local (Smith-Waterman, BLOSUM62) alignment of that peptide
within each protein sequence. This recovers a match even where the real
protein differs from every one of the primer's encoded possibilities --
e.g. a substitution the primer wasn't designed to tolerate, a slightly
divergent homolog, or a sequencing error -- at the cost of being less
exact than a true regex match: a local alignment always returns its best
attempt, even against a completely unrelated sequence, so the reported
score/identity for each hit is worth checking rather than trusting blindly
(``--min_identity`` can filter low-confidence hits automatically).

For each sequence, the "amplicon" is the fixed-length window of
``--amplicon_size`` amino acids (33 by default, matching the amplicon size
used throughout the rest of the ``amplicons`` package) starting at that
anchor -- not necessarily where the reported local alignment happens to
start. A local alignment trims off leading residues whenever including
them would lower the score (e.g. a mismatched first residue), even though
the target does have *some* residue there; anchoring to the trimmed start
would silently drop that real upstream sequence from the amplicon. Instead
the anchor is walked back by however many leading consensus-peptide
residues got trimmed, so the amplicon covers the same span the primer
would regardless of whether its very first residue happened to align. See
``find_best_alignment_amplicon`` for the details and its assumptions.
"""

from argparse import ArgumentParser, Namespace
from dataclasses import dataclass
from typing import Optional

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from poliwag.amplicons.degenerate_primer import translate_codons

# The degenerate forward primer this script was written for. Pass a
# different primer with -p/--primer to reuse it for another assay.
DEFAULT_PRIMER = "GCSTACSYSATSTACACSTCSGG"

# Amplicon length (in amino acids) extracted starting at each alignment,
# matching the amplicon size used elsewhere in the amplicons package.
DEFAULT_AMPLICON_SIZE = 33

# EMBOSS water/needle-style affine gap penalties for protein alignment.
DEFAULT_OPEN_GAP_SCORE = -10.0
DEFAULT_EXTEND_GAP_SCORE = -0.5


def parse_arguments() -> Namespace:
    parser = ArgumentParser(
        description="Locally align a consensus peptide translated from a degenerate DNA primer "
        "against a FASTA file of protein sequences, and extract the fixed-length amplicon "
        "starting where each best alignment begins."
    )
    parser.add_argument("-f", "--fasta", type=str, required=True,
                        help="Path to FASTA file containing the protein sequences to search")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Path to output FASTA file of extracted amplicons")
    parser.add_argument("-p", "--primer", type=str, default=DEFAULT_PRIMER,
                        help=f"Degenerate DNA primer, IUPAC ambiguity codes allowed "
                        f"(default: {DEFAULT_PRIMER})")
    parser.add_argument("-a", "--amplicon_size", type=int, default=DEFAULT_AMPLICON_SIZE,
                        help=f"Length in amino acids of the amplicon extracted starting at each "
                        f"alignment (default: {DEFAULT_AMPLICON_SIZE})")
    parser.add_argument("--min_identity", type=float, default=0.0,
                        help="Skip alignments below this percent identity over the aligned region "
                        "(0-100, default: 0, i.e. report every alignment and filter afterwards)")
    parser.add_argument("--open_gap_score", type=float, default=DEFAULT_OPEN_GAP_SCORE,
                        help=f"Gap opening penalty for the local alignment (default: {DEFAULT_OPEN_GAP_SCORE})")
    parser.add_argument("--extend_gap_score", type=float, default=DEFAULT_EXTEND_GAP_SCORE,
                        help=f"Gap extension penalty for the local alignment (default: {DEFAULT_EXTEND_GAP_SCORE})")
    args = parser.parse_args()
    return args


def primer_consensus_peptide(primer: str) -> str:
    """Translate a degenerate DNA primer to a single representative peptide.

    For each of the primer's codons, picks the single most common amino
    acid it encodes across all of its concrete (non-degenerate) codons.
    Ties are broken alphabetically, so this is one of possibly several
    equally "most common" peptides the primer could encode -- a single
    representative to align with, not a claim about the true original
    amino acid at each position.

    :param primer: degenerate DNA primer sequence
    :return: consensus peptide translation of the primer's complete codons
    """
    translations, _leftover = translate_codons(primer)
    return "".join(t.consensus for t in translations)


def build_aligner(open_gap_score: float, extend_gap_score: float) -> PairwiseAligner:
    """Configure a local (Smith-Waterman) protein aligner.

    :param open_gap_score: gap opening penalty (negative)
    :param extend_gap_score: gap extension penalty (negative)
    :return: a ``PairwiseAligner`` in local mode, scored with BLOSUM62
    """
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score
    return aligner


@dataclass
class AlignmentAmplicon:
    """The best local alignment of the consensus peptide to one protein
    sequence, and its extracted amplicon."""

    record_id: str
    anchor_start: int  # 0-based start of the amplicon in the target (see find_best_alignment_amplicon)
    aligned_target_start: int  # 0-based start of the *aligned region itself* in the target (>= anchor_start)
    aligned_target_end: int  # 0-based, exclusive end of the aligned region in the target
    score: float
    identity: float  # percent identity over the aligned region
    amplicon: str
    trimmed_upstream: int  # leading consensus-peptide residues the local alignment trimmed off


def find_best_alignment_amplicon(
    record: SeqRecord, consensus_peptide: str, aligner: PairwiseAligner, amplicon_size: int
) -> Optional[AlignmentAmplicon]:
    """Locally align ``consensus_peptide`` to ``record`` and extract the
    fixed-length amplicon anchored to where the *full* consensus peptide
    would start in the target -- not necessarily where the reported
    alignment happens to start.

    A local (Smith-Waterman) alignment trims off leading/trailing residues
    whenever including them would lower the score -- e.g. if the target's
    first residue is a mismatch against the consensus peptide's first
    residue, the alignment simply won't include it, and will report
    starting one residue later even though the target does have a
    (mismatched) residue there. Anchoring the amplicon to that trimmed
    start would silently drop real upstream sequence. Instead, the anchor
    is walked back by however many leading consensus-peptide residues were
    trimmed (``trimmed_upstream``), so the amplicon covers the same span
    the primer would regardless of whether its very first residue happened
    to align.

    This assumes a 1:1 correspondence over that trimmed leading stretch
    (i.e. no gap hiding in it), which holds for a short, gap-free
    consensus peptide like this one aligned with standard affine gap
    penalties -- a gap this close to the alignment's edge is never
    favourable over simply trimming.

    ``aligner.align`` is only asked for one alignment (via ``next(iter(...))``)
    rather than all optimal alignments, since a local alignment can have many
    equally-scoring tracebacks and enumerating them all is unnecessary here
    -- any one of the best-scoring alignments gives the same anchor for a
    short, gap-free consensus peptide like this one.

    :param record: protein sequence to align against
    :param consensus_peptide: representative peptide translated from the primer
    :param aligner: configured local ``PairwiseAligner``
    :param amplicon_size: length in amino acids of the amplicon to extract,
        starting at the anchor described above
    :return: the best-scoring alignment's amplicon, or ``None`` if the
        target doesn't have ``amplicon_size`` residues remaining from the
        anchor (this includes the target being too short upstream to walk
        the anchor back by ``trimmed_upstream``)
    :raises ValueError: if ``record`` contains a character outside the
        substitution matrix's alphabet (e.g. a gap "-", a lowercase
        letter, a digit, or "U"/selenocysteine) -- callers should catch
        this per-record rather than let it abort the whole run
    """
    target_sequence = str(record.seq)
    alignment = next(iter(aligner.align(consensus_peptide, target_sequence)))

    trimmed_upstream = alignment.aligned[0][0][0]  # leading consensus residues left out of the alignment
    aligned_target_start = alignment.aligned[1][0][0]
    aligned_target_end = alignment.aligned[1][-1][-1]
    anchor_start = aligned_target_start - trimmed_upstream

    if anchor_start < 0:
        return None

    amplicon_end = anchor_start + amplicon_size
    if amplicon_end > len(target_sequence):
        return None

    # counts.aligned isn't present in every Biopython version (it's a newer
    # addition to AlignmentCounts), so the aligned-column count is derived
    # from identities + mismatches instead, which has been available for
    # much longer and is equivalent for an ungapped consensus peptide like
    # this one.
    counts = alignment.counts()
    aligned_columns = counts.identities + counts.mismatches
    identity = 100 * counts.identities / aligned_columns if aligned_columns else 0.0

    return AlignmentAmplicon(
        record_id=record.id,
        anchor_start=anchor_start,
        aligned_target_start=aligned_target_start,
        aligned_target_end=aligned_target_end,
        score=alignment.score,
        identity=identity,
        amplicon=target_sequence[anchor_start:amplicon_end],
        trimmed_upstream=trimmed_upstream,
    )


def main() -> None:
    args = parse_arguments()

    consensus_peptide = primer_consensus_peptide(args.primer)
    translations, leftover = translate_codons(args.primer)
    print(f"Primer: {args.primer}")
    print(f"Consensus peptide: {consensus_peptide}")
    for i, translation in enumerate(translations, start=1):
        if translation.is_tie:
            print(f"  codon {i} ({translation.codon}): tied between "
                  f"{''.join(translation.possible)}, resolved to "
                  f"'{translation.consensus}' (alphabetically first)")
    if leftover:
        print(f"Note: last {len(leftover)} base(s) of the primer ({leftover}) don't complete a "
              f"codon and were ignored.")

    aligner = build_aligner(args.open_gap_score, args.extend_gap_score)

    amplicons_written = 0
    records_skipped_short = 0
    records_skipped_identity = 0
    records_skipped_invalid = 0

    with open(args.output, "w") as out:
        for record in tqdm(SeqIO.parse(args.fasta, "fasta"), desc="Aligning sequences", unit="sequence"):
            try:
                hit = find_best_alignment_amplicon(record, consensus_peptide, aligner, args.amplicon_size)
            except ValueError as error:
                # Raised by PairwiseAligner when the sequence contains a
                # character outside BLOSUM62's alphabet (a gap "-", a
                # lowercase letter, a digit, "U"/selenocysteine, etc). Warn
                # and move on rather than letting one bad sequence crash
                # the whole run.
                print(f"Warning: skipping {record.id}: {error}")
                records_skipped_invalid += 1
                continue
            if hit is None:
                records_skipped_short += 1
                continue
            if hit.identity < args.min_identity:
                records_skipped_identity += 1
                continue

            header = (f"{hit.record_id}|primer_alignment_{hit.anchor_start}-"
                      f"{hit.anchor_start + args.amplicon_size}|score={hit.score:.1f}"
                      f"|identity={hit.identity:.1f}%")
            if hit.trimmed_upstream:
                header += f"|trimmed_upstream={hit.trimmed_upstream}"
            out.write(f">{header}\n{hit.amplicon}\n")
            amplicons_written += 1

    print(f"Wrote {amplicons_written} amplicon(s).")
    print(f"{records_skipped_short} sequence(s) skipped: not enough sequence for a full "
          f"{args.amplicon_size}-residue amplicon, once anchored to where the full consensus "
          f"peptide would start.")
    print(f"{records_skipped_identity} sequence(s) skipped: below {args.min_identity}% identity.")
    print(f"{records_skipped_invalid} sequence(s) skipped: contained a character the aligner "
          f"couldn't score (see warnings above).")


if __name__ == "__main__":
    main()
