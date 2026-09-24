"""
Downstream analysis of a NATU cluster: link a query polymer's cluster to its real
genome occurrences, split those occurrences into "matches" and "variants" of the
query polymer, and use AMP-binding domain amplicon similarity to characterize how
distinct those groups are from each other and from other clusters.

Pipeline assumed (matching the rest of the poliwag/NATU workflow):

1. A FASTA file of *unique* polymers was built (e.g. ``polymers_to_fasta.py``),
   with arbitrary headers (e.g. "polymer_1") and one pipe-joined monomer sequence
   per record.
2. That FASTA was clustered with ``natu cluster``, producing a ``clusters.tsv``
   (columns: header, cluster_id) via ``natu.network.write_clusters``.
3. Real genome occurrences of those unique polymers live in a polymer metadata
   file (e.g. ``natu_all_polymers.txt``, columns: Accession, Region, Candidate,
   Polymer, Domains), readable with
   ``poliwag.antismash_processing.get_polymers.read_polymers``. The same unique
   polymer sequence can appear in this file multiple times (multiple genomes/
   regions), which is why one unique polymer maps to *several* occurrences. The
   "Domains" column (see ``get_polymers.write_polymers``) records, for each
   monomer, the (locus_tag, domain_index) of the AMP-binding domain that
   produced it -- this is what makes position-paired amplicon comparison
   (below) possible; a metadata file written before that column existed can
   still be read, but its occurrences simply can't take part in any
   comparison (there's no way to know which amplicon is which domain).
4. AMP-binding domain amplicons were extracted with
   ``poliwag.amplicons.amp_binding_domains_to_fasta`` and then reduced to fixed
   windows with ``poliwag.amplicons.in_silico_amplicon_sequencing``, giving a
   FASTA whose headers start with
   ``<accession>|region_<region>|candidate_<candidates>|<locus_tag>|<domain_index>``
   -- the same accession/region/candidate/locus_tag/domain_index identifiers
   used in the polymer metadata file, which is what lets a specific polymer
   *position* (not just a whole occurrence) be resolved to its one amplicon.

Given a query polymer (as a pipe-joined monomer string) and these files, this
script:

1. Finds the NATU cluster containing the query polymer, and all unique polymer
   members of that cluster.
2. Resolves every unique member to all of its real genome occurrences in the
   polymer metadata file (exact sequence match; a unique polymer can have many
   occurrences).
3. Splits those unique members into "matches" (a member whose own polymer
   sequence exactly equals the query polymer's sequence, or contains it as a
   contiguous run of monomers) and "variants" (everything else in the
   cluster).
4. For every comparison -- within matches, matches vs. variants, and matches
   vs. one representative occurrence per *other* cluster ("non-matches") --
   amplicons are compared only at *corresponding polymer positions*, never as
   a pooled bag of "every domain against every other domain":
   - For a match, the query occurs verbatim as a contiguous run, so its
     query-position -> match-position correspondence is a plain offset (e.g.
     query position 0 is match position ``offset``).
   - For a variant or a non-match representative, the two polymers can
     differ by substitutions/insertions/deletions, so the correspondence is
     found by globally aligning the two monomer sequences (treating each
     monomer name as one token, via ``Bio.Align.PairwiseAligner``) and
     keeping only the aligned (non-gap) column pairs. An inserted or deleted
     monomer on either side simply contributes no comparison at that
     position, rather than being forced into one. If ``--search_tsv`` is
     given (a precomputed ``natu search`` output TSV), an already-computed
     alignment is reused from it instead whenever that exact query/subject
     pair is present there; see ``load_search_alignments`` for what that
     requires.
   This is what lets "the Trp-A domain" only ever get compared against
   "the Trp-A-equivalent domain" in the other polymer, never against an
   unrelated Asn-A or Thr-A domain from the same occurrence.
5. Reports both a single blended average per comparison (for a quick read)
   and a per-polymer-position breakdown against the query (so a position
   that behaves very differently from the rest -- more variable, or better
   conserved -- doesn't get hidden inside a blended number).
6. Prints a suggested percent-identity range for calling something a
   "predicted variant" of the query polymer, derived from where the observed
   variant similarity sits between the match and non-match baselines. This is
   a heuristic based on a single cluster's data -- treat it as a starting
   point, not a validated cutoff, until checked against more clusters.
7. Optionally (``--verified_amplicon_fasta``), also compares a verified set
   of amplicons for this cluster -- one sequence per query polymer position,
   in polymer order, not tied to any single genome occurrence -- against the
   *in silico* (automated) amplicons of the match occurrences, the variant
   occurrences, and the non-match representatives, the same position-paired
   way as step 4. This asks a different question from step 4's
   occurrence-vs-occurrence comparisons: not "how similar are these groups
   of automated calls to each other", but "how well does the automated
   in-silico amplicon extraction agree with a trusted reference, for each
   group" -- e.g. it should score highest for matches, lower for variants,
   and lowest for non-matches, sanity-checking the in-silico extraction
   itself rather than only the polymers it was run on.
"""

import csv
import random
from argparse import ArgumentParser, Namespace
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path
from statistics import mean
from typing import Iterable, Optional

import numpy as np
from Bio.Align import PairwiseAligner

from poliwag.antismash_processing.get_polymers import read_polymers
from poliwag.antismash_processing.polymer import PredictedPolymer
from poliwag.utils import parse_fasta_file


# ---------------------------------------------------------------------------
# Amplicon percent identity
# ---------------------------------------------------------------------------

_IDENTITY_ALIGNER = PairwiseAligner()
_IDENTITY_ALIGNER.mode = "global"
_IDENTITY_ALIGNER.match_score = 1.0
_IDENTITY_ALIGNER.mismatch_score = 0.0
_IDENTITY_ALIGNER.open_gap_score = 0.0
_IDENTITY_ALIGNER.extend_gap_score = 0.0


def percent_identity(seq1: str, seq2: str) -> float:
    """Percent identity between two amplicon (protein) sequences.

    Amplicons are normally the same fixed length (they're windows extracted at
    a common HMM alignment column), so the common case is a plain positional
    comparison. If the two sequences differ in length (e.g. a truncated
    amplicon near a sequence end), falls back to a simple global alignment
    (match=1, mismatch=0, no gap penalty) and divides by the longer sequence's
    length, so a partial/truncated amplicon doesn't get an inflated score.

    :param seq1: first amplicon sequence
    :param seq2: second amplicon sequence
    :return: percent identity, 0-100
    """
    if len(seq1) == len(seq2):
        if not seq1:
            return 0.0
        matches = sum(a == b for a, b in zip(seq1, seq2))
        return 100.0 * matches / len(seq1)

    alignment_score = _IDENTITY_ALIGNER.score(seq1, seq2)
    return 100.0 * alignment_score / max(len(seq1), len(seq2))


def average_paired_identity(pairs: list[tuple[str, str]]) -> Optional[float]:
    """Average percent identity over an explicit list of already
    position-paired amplicon sequences (see the module docstring -- this is
    deliberately *not* an all-vs-all pool; every pair here is expected to
    come from corresponding polymer positions).

    :param pairs: list of (sequence_a, sequence_b) amplicon pairs
    :return: average percent identity, or None if ``pairs`` is empty
    """
    if not pairs:
        return None
    return mean(percent_identity(a, b) for a, b in pairs)


# ---------------------------------------------------------------------------
# Polymer sequence matching
# ---------------------------------------------------------------------------

def monomers_of(polymer_string: str, separator: str) -> list[str]:
    return polymer_string.split(separator)


def contains_subsequence(haystack: list[str], needle: list[str]) -> Optional[int]:
    """The start index of ``needle`` as a contiguous run within ``haystack``,
    or None if it doesn't appear as one. The first such run is used if there
    happens to be more than one (rare for real NRPS assembly lines)."""
    if len(needle) > len(haystack):
        return None
    for i in range(len(haystack) - len(needle) + 1):
        if haystack[i:i + len(needle)] == needle:
            return i
    return None


def match_offset(candidate_monomers: list[str], query_monomers: list[str]) -> Optional[int]:
    """Whether ``candidate_monomers`` is a "match" to ``query_monomers`` --
    identical to it, or containing it as a contiguous run -- and if so, the
    index in ``candidate_monomers`` where the query starts (0 for an exact
    match). Returns None if it isn't a match at all.

    This offset is what lets every position in a match occurrence be
    resolved back to "the same slot" in the query without needing any
    alignment: query position i is always candidate position ``offset + i``.
    """
    if candidate_monomers == query_monomers:
        return 0
    return contains_subsequence(candidate_monomers, query_monomers)


def is_match(candidate_monomers: list[str], query_monomers: list[str]) -> bool:
    """Whether ``candidate_monomers`` is a "match" to ``query_monomers`` (see
    ``match_offset``); kept as a readable yes/no wrapper where the offset
    itself isn't needed."""
    return match_offset(candidate_monomers, query_monomers) is not None


# ---------------------------------------------------------------------------
# Monomer-level alignment (for variants / non-matches, where there's no
# simple offset relating the two polymers)
# ---------------------------------------------------------------------------

_MONOMER_ALIGNER = PairwiseAligner()
_MONOMER_ALIGNER.mode = "global"
_MONOMER_ALIGNER.match_score = 1.0
_MONOMER_ALIGNER.mismatch_score = 0.0
_MONOMER_ALIGNER.open_gap_score = -1.0
_MONOMER_ALIGNER.extend_gap_score = -0.5


def _encode_monomers(query_monomers: list[str], other_monomers: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """Encode two monomer-name sequences as int32 arrays over a shared alphabet
    (the union of both sequences' distinct monomer names), for feeding to
    Biopython's ``PairwiseAligner`` instead of the raw string lists.

    This sidesteps a real Biopython compatibility problem: whether
    ``aligner.align()`` accepts a plain ``list[str]`` of arbitrary (not
    necessarily single-character) tokens at all -- and if so, whether it can
    infer an alphabet for one on its own -- has changed across Biopython
    versions. Confirmed empirically: passing ``list[str]`` monomer sequences
    directly raises ``ValueError: alphabet is None; cannot interpret
    sequence`` on Biopython 1.79/1.81, and ``TypeError: sequence has
    unexpected type list`` on 1.83; both reproduce with nothing more unusual
    than ordinary multi-word monomer names (e.g. "aspartic acid"), not just
    some degenerate input. It happens to work on 1.88, but relying on that
    is exactly the kind of thing that breaks again, silently, on someone
    else's (or a future) environment. Pre-converting to int arrays with an
    explicit alphabet -- the same way ``natu.pairwise.Converter`` does for
    NATU's own alignments -- works the same way across all of these
    versions; confirmed on 1.79 specifically, which is where this was first
    hit (as ``ValueError: alphabet is None; cannot interpret sequence``).

    The two returned arrays are only meaningful compared against each other
    (same alphabet, built from their own union) -- not against arrays
    encoded from any other pair.

    :param query_monomers: first monomer sequence
    :param other_monomers: second monomer sequence
    :return: ``(query_codes, other_codes)``, both int32 numpy arrays
    """
    alphabet = sorted(set(query_monomers) | set(other_monomers))
    index = {monomer: code for code, monomer in enumerate(alphabet)}
    query_codes = np.array([index[m] for m in query_monomers], dtype=np.int32)
    other_codes = np.array([index[m] for m in other_monomers], dtype=np.int32)
    return query_codes, other_codes


def align_monomer_positions(query_monomers: list[str], other_monomers: list[str]) -> list[tuple[int, int]]:
    """Globally align two polymers at the monomer level (each monomer name is
    one token, not one character -- encoded to an int array first, see
    ``_encode_monomers``, so Biopython's ``PairwiseAligner`` treats each
    monomer name as one indivisible symbol) and return the list of
    (query_index, other_index) pairs for every aligned column, in polymer
    order.

    An aligned pair can be a match *or* a mismatch (a substituted monomer) --
    both are kept, since the point is comparing amplicon identity at
    "the same slot", regardless of whether the monomer itself changed. Only
    genuinely unpaired (inserted/deleted) monomers are left out.
    """
    if not query_monomers or not other_monomers:
        return []
    query_codes, other_codes = _encode_monomers(query_monomers, other_monomers)
    alignments = _MONOMER_ALIGNER.align(query_codes, other_codes)
    if not alignments:
        return []
    query_blocks, other_blocks = alignments[0].aligned
    pairs: list[tuple[int, int]] = []
    for (q_start, q_end), (o_start, o_end) in zip(query_blocks, other_blocks):
        for step in range(q_end - q_start):
            pairs.append((q_start + step, o_start + step))
    return pairs


# ---------------------------------------------------------------------------
# Precomputed alignments from `natu search` output (optional)
# ---------------------------------------------------------------------------

# NATU's own gap sentinel in a `natu search`/`natu align` TSV's aligned_* columns
# (natu.constants.GAP_REPR) -- duplicated here, not imported, so this script stays
# decoupled from the natu package, the same way load_clusters() reads clusters.tsv
# as plain TSV instead of importing natu.network.
_NATU_GAP_REPR = "GAP"
# Monomer join separator NATU itself always uses in its own FASTA/alignment output
# (see the natu-project-reference: "one pipe-joined monomer sequence per record") --
# fixed, independent of this script's own --separator argument.
_NATU_MONOMER_SEPARATOR = "|"


def _decode_search_alignment_row(
    aligned_query: str, aligned_subject: str
) -> tuple[tuple[str, ...], tuple[str, ...], list[tuple[int, int]]]:
    """Decode one `natu search` output row's ``aligned_q``/``s_aligned_s`` columns
    (equal-length, pipe-joined, gap-padded monomer sequences) into the ungapped
    query/subject monomer tuples plus the (query_index, subject_index) pairs for
    every aligned (non-gap-non-gap) column -- the same shape
    ``align_monomer_positions`` returns for this pair, read off disk instead of
    recomputed.
    """
    query_tokens = aligned_query.split(_NATU_MONOMER_SEPARATOR)
    subject_tokens = aligned_subject.split(_NATU_MONOMER_SEPARATOR)

    query_monomers: list[str] = []
    subject_monomers: list[str] = []
    position_pairs: list[tuple[int, int]] = []

    for query_token, subject_token in zip(query_tokens, subject_tokens):
        query_is_gap = query_token == _NATU_GAP_REPR
        subject_is_gap = subject_token == _NATU_GAP_REPR
        if not query_is_gap and not subject_is_gap:
            position_pairs.append((len(query_monomers), len(subject_monomers)))
        if not query_is_gap:
            query_monomers.append(query_token)
        if not subject_is_gap:
            subject_monomers.append(subject_token)

    return tuple(query_monomers), tuple(subject_monomers), position_pairs


def load_search_alignments(
    search_tsv: str,
) -> dict[tuple[tuple[str, ...], tuple[str, ...]], list[tuple[int, int]]]:
    """Load a `natu search` output TSV (columns: query, subject, bitscore, aligned_q,
    s_aligned_s -- see ``natu.cli``'s search writer) into a lookup from
    ``(query_monomers, other_monomers)`` (as tuples of monomer names, so lookups
    don't depend on this script's own ``--separator`` or on whatever FASTA headers
    NATU happened to use) to the (query_index, other_index) position pairs for every
    aligned column -- the same shape ``align_monomer_positions`` returns, so callers
    can use whichever is available without caring which one it came from.

    For a query/subject pair to be found here, it must have been an actual row in
    this file -- i.e. the `natu search` run that produced it must have scored that
    exact pair above its own ``--threshold``. Because NATU's clustering is
    single-linkage (see the natu-project-reference skill), a genuine cluster member
    is not guaranteed to score above any nonzero threshold directly against the
    query, so for full coverage of every cluster member -- and of whichever
    occurrence the non-match baseline happens to randomly sample -- generate this
    file with ``natu search --threshold 0`` against the same ``--cluster_fasta``
    used here, ideally with the same substitution matrix/config as the ``natu
    cluster`` run that produced ``--clusters``. Any pair missing from this file is
    simply absent from the returned dict; callers fall back to computing it fresh
    (see ``resolve_variant_position_pairs``).

    If the same (query, subject) monomer pair appears in more than one row, the
    first row wins and later duplicates are ignored.

    :param search_tsv: path to the `natu search` output TSV
    """
    alignments: dict[tuple[tuple[str, ...], tuple[str, ...]], list[tuple[int, int]]] = {}
    with open(search_tsv, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            query_monomers, subject_monomers, position_pairs = _decode_search_alignment_row(
                row["aligned_q"], row["s_aligned_s"]
            )
            alignments.setdefault((query_monomers, subject_monomers), position_pairs)
    return alignments


def resolve_variant_position_pairs(
    query_monomers: list[str],
    other_monomers: list[str],
    search_alignments: Optional[dict[tuple[tuple[str, ...], tuple[str, ...]], list[tuple[int, int]]]],
) -> tuple[list[tuple[int, int]], bool]:
    """Position pairs for ``other_monomers`` against ``query_monomers``: reused from
    ``search_alignments`` (see ``load_search_alignments``) when it's given and has
    this exact pair, otherwise computed fresh with ``align_monomer_positions``.

    :return: ``(position_pairs, from_search)`` -- ``from_search`` is True if the
        pairs came from ``search_alignments`` rather than a fresh alignment.
    """
    if search_alignments is not None:
        cached = search_alignments.get((tuple(query_monomers), tuple(other_monomers)))
        if cached is not None:
            return cached, True
    return align_monomer_positions(query_monomers, other_monomers), False


# ---------------------------------------------------------------------------
# Amplicon lookup
# ---------------------------------------------------------------------------

@dataclass
class AmpliconRecord:
    accession: str
    region: int
    candidate_clusters: set[int]
    locus_tag: str
    domain_index: str
    sequence: str


def _parse_candidate_field(field_str: str) -> set[int]:
    """Parse a 'candidate_<...>' fasta-id field into a set of ints.

    Handles the plain case ("candidate_19" -> {19}), the interleaved/
    overlapping case ("candidate_23-24" -> {23, 24}), and the no-candidate
    case ("candidate_NA" -> set()).
    """
    value = field_str.split("_", 1)[1] if "_" in field_str else field_str
    if value == "NA":
        return set()
    return {int(v) for v in value.split("-")}


def load_amplicons(amplicon_fasta: str, separator: str = "|") -> list[AmpliconRecord]:
    """Load amplicon sequences, parsing the accession/region/candidate/locus_tag
    identifiers back out of each header (as written by
    ``poliwag.amplicons.amp_binding_domains_to_fasta`` and carried through
    ``poliwag.amplicons.in_silico_amplicon_sequencing``).

    :param amplicon_fasta: path to the amplicon FASTA file
    :param separator: separator used between id fields in the header (not the
        monomer separator)
    """
    records: list[AmpliconRecord] = []
    fasta_dict = parse_fasta_file(amplicon_fasta)
    for header, sequence in fasta_dict.items():
        fields = header.split(separator)
        if len(fields) < 5:
            print(f"Warning: could not parse amplicon header {header!r}, skipping")
            continue
        accession, region_field, candidate_field, locus_tag, domain_index = fields[:5]
        try:
            region = int(region_field.split("_", 1)[1])
            candidate_clusters = _parse_candidate_field(candidate_field)
        except (IndexError, ValueError):
            print(f"Warning: could not parse amplicon header {header!r}, skipping")
            continue

        records.append(AmpliconRecord(
            accession=accession,
            region=region,
            candidate_clusters=candidate_clusters,
            locus_tag=locus_tag,
            domain_index=domain_index,
            sequence=sequence,
        ))
    return records


def load_verified_amplicons(verified_amplicon_fasta: str, expected_length: int) -> list[str]:
    """Load a FASTA of verified AMP-binding domain amplicons for the query's
    own cluster: one record per polymer position, in file order (position 0
    first) -- not tied to any accession/region/locus_tag/domain_index the way
    ``load_amplicons`` is, since this represents a single trusted reference
    for the query polymer itself, not a pool of per-occurrence amplicons.
    Headers are not parsed or otherwise interpreted; only each record's order
    and its sequence matter.

    :param verified_amplicon_fasta: path to the verified amplicon FASTA
    :param expected_length: the query polymer's monomer count -- checked
        against the number of records here, since a mismatch would silently
        misalign every downstream position-paired comparison against it
        rather than raising anywhere obvious
    :raises ValueError: if the record count doesn't match ``expected_length``
    """
    fasta_dict = parse_fasta_file(verified_amplicon_fasta)
    verified_amplicons = list(fasta_dict.values())
    if len(verified_amplicons) != expected_length:
        raise ValueError(
            f"--verified_amplicon_fasta has {len(verified_amplicons)} record(s), but the "
            f"query polymer has {expected_length} monomer(s) -- these must correspond "
            "1:1, in polymer order, for a position-paired comparison against it to mean "
            "anything."
        )
    return verified_amplicons


def index_amplicons_by_domain(amplicons: list[AmpliconRecord]) -> dict[tuple[str, int, str, str], AmpliconRecord]:
    """Index amplicons by (accession, region, locus_tag, domain_index) -- the
    exact key stored in ``PredictedPolymer.domain_ids`` -- for O(1) lookup of
    the one amplicon that corresponds to a specific polymer position."""
    index: dict[tuple[str, int, str, str], AmpliconRecord] = {}
    for amplicon in amplicons:
        key = (amplicon.accession, amplicon.region, amplicon.locus_tag, amplicon.domain_index)
        index[key] = amplicon
    return index


def ordered_amplicons_for_occurrence(
    occurrence: PredictedPolymer,
    amplicon_by_domain: dict[tuple[str, int, str, str], AmpliconRecord],
) -> Optional[list[Optional[AmpliconRecord]]]:
    """The amplicon at each polymer position of ``occurrence``, in the same
    order as ``occurrence.polymer`` -- i.e. entry i of the returned list is
    the AMP-binding domain amplicon that encodes ``occurrence.polymer[i]``,
    or None if that specific position has no corresponding amplicon (e.g. a
    PKS_AT/CAL-derived monomer, or a missing PARAS/hmmscan result).

    Returns None (not a list) if ``occurrence.domain_ids`` wasn't recorded at
    all -- e.g. the polymer metadata file was written before the "Domains"
    column existed. There's no way to do position-paired comparison for such
    an occurrence; callers should skip it (and it's worth re-running
    ``get_polymers`` on a version of antiSMASH/poliwag that writes it).
    """
    if occurrence.domain_ids is None:
        return None
    ordered: list[Optional[AmpliconRecord]] = []
    for domain_id in occurrence.domain_ids:
        if domain_id is None:
            ordered.append(None)
            continue
        locus_tag, domain_index = domain_id
        key = (occurrence.accession, occurrence.region, locus_tag, domain_index)
        ordered.append(amplicon_by_domain.get(key))
    return ordered


def amplicon_at(ordered: Optional[list[Optional[AmpliconRecord]]], index: int) -> Optional[AmpliconRecord]:
    """Safe positional lookup into an ``ordered_amplicons_for_occurrence``
    result: None if there's no provenance at all, or the index is out of
    range, or that specific position has no amplicon."""
    if ordered is None or not (0 <= index < len(ordered)):
        return None
    return ordered[index]


def accumulate_by_query_position(
    by_position: dict[int, list[tuple[str, str]]],
    ordered_a: Optional[list[Optional[AmpliconRecord]]],
    ordered_b: Optional[list[Optional[AmpliconRecord]]],
    triples: Iterable[tuple[int, int, int]],
) -> None:
    """For each ``(query_position, index_in_a, index_in_b)`` triple, resolve
    the amplicon at ``index_in_a`` in ``ordered_a`` and at ``index_in_b`` in
    ``ordered_b``, and -- if both sides actually have one -- append the pair
    to ``by_position[query_position]`` (mutated in place, so callers can
    accumulate across many occurrence pairs into one shared dict).

    ``query_position`` is tracked separately from the two lookup indices
    because they can differ from it and from each other: for a match
    occurrence the lookup index is ``offset + query_position`` (a plain
    shift), and for a variant/non-match occurrence it's whatever
    ``align_monomer_positions`` mapped that query position to.
    """
    for query_position, index_a, index_b in triples:
        amp_a = amplicon_at(ordered_a, index_a)
        amp_b = amplicon_at(ordered_b, index_b)
        if amp_a is not None and amp_b is not None:
            by_position.setdefault(query_position, []).append((amp_a.sequence, amp_b.sequence))


def accumulate_verified_by_position(
    by_position: dict[int, list[tuple[str, str]]],
    verified_amplicons: list[str],
    ordered_other: Optional[list[Optional[AmpliconRecord]]],
    pairs: Iterable[tuple[int, int]],
) -> None:
    """Like ``accumulate_by_query_position``, but one side is fixed: the
    verified amplicon at each query polymer position (see
    ``load_verified_amplicons``), rather than a second occurrence's amplicon.

    ``pairs`` is ``(query_index, other_index)`` -- ``other_index`` is
    resolved against ``ordered_other`` the same way as elsewhere; the
    verified side is indexed directly by ``query_index`` since
    ``verified_amplicons`` is already in query-polymer order and has no
    "occurrence has no amplicon at all" case the way a real genome
    occurrence can.
    """
    for query_index, other_index in pairs:
        amp_other = amplicon_at(ordered_other, other_index)
        if amp_other is not None:
            by_position.setdefault(query_index, []).append(
                (verified_amplicons[query_index], amp_other.sequence)
            )


def split_position_pairs_by_substrate(
    query_monomers: list[str],
    other_monomers: list[str],
    position_pairs: Iterable[tuple[int, int]],
) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    """Split ``(query_index, other_index)`` position pairs -- e.g. a variant
    group's monomer-level alignment, see ``_VariantGroup.position_pairs`` --
    into "substrate match" pairs, where the other polymer's monomer at
    ``other_index`` is the *same* substrate as the query's monomer at
    ``query_index`` (this specific position wasn't substituted, even though
    the polymer as a whole differs elsewhere), and "substrate mismatch"
    pairs, where it's a genuinely different monomer (an actual substitution
    at this position).

    :param query_monomers: the query polymer's monomers, indexed by
        ``query_index``
    :param other_monomers: the other polymer's monomers (e.g. a
        ``_VariantGroup.polymer``), indexed by ``other_index``
    :param position_pairs: the ``(query_index, other_index)`` pairs to split
    :return: ``(substrate_match_pairs, substrate_mismatch_pairs)``, each a
        list of the same kind of ``(query_index, other_index)`` pairs,
        partitioned by whether the two monomers are equal
    """
    substrate_match_pairs: list[tuple[int, int]] = []
    substrate_mismatch_pairs: list[tuple[int, int]] = []
    for query_index, other_index in position_pairs:
        if query_monomers[query_index] == other_monomers[other_index]:
            substrate_match_pairs.append((query_index, other_index))
        else:
            substrate_mismatch_pairs.append((query_index, other_index))
    return substrate_match_pairs, substrate_mismatch_pairs


# ---------------------------------------------------------------------------
# NATU cluster loading
# ---------------------------------------------------------------------------

def load_clusters(clusters_tsv: str) -> dict[str, int]:
    """Load a NATU ``clusters.tsv`` (as written by ``natu.network.write_clusters``)
    into a header -> cluster_id mapping."""
    header_to_cluster: dict[str, int] = {}
    with open(clusters_tsv, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            header_to_cluster[row["header"]] = int(row["cluster_id"])
    return header_to_cluster


@dataclass
class PositionResult:
    query_index: int
    query_monomer: str
    avg_identity: Optional[float]
    n_pairs: int


@dataclass
class ClusterAnalysis:
    query_polymer: str
    query_monomers: list[str]
    cluster_id: int
    cluster_members: list[str] = field(default_factory=list)  # unique polymer strings
    match_occurrences: list[PredictedPolymer] = field(default_factory=list)
    variant_occurrences: list[PredictedPolymer] = field(default_factory=list)
    occurrences_without_provenance: int = 0  # match/variant occurrences skipped: no Domains data
    avg_within_matches: Optional[float] = None
    avg_matches_variants: Optional[float] = None
    avg_matches_nonmatches: Optional[float] = None
    n_within_matches_pairs: int = 0
    n_matches_variants_pairs: int = 0
    n_matches_nonmatches_pairs: int = 0
    within_matches_by_position: list[PositionResult] = field(default_factory=list)
    matches_variants_by_position: list[PositionResult] = field(default_factory=list)
    matches_nonmatches_by_position: list[PositionResult] = field(default_factory=list)
    non_match_cluster_ids: list[int] = field(default_factory=list)
    skipped_non_match_clusters: list[int] = field(default_factory=list)
    used_search_tsv: bool = False  # whether a --search_tsv was supplied at all
    n_position_pairs_from_search_tsv: int = 0  # variant/non-match alignments reused from it
    n_position_pairs_computed: int = 0  # ...and how many fell back to a fresh alignment
    used_verified_amplicons: bool = False  # whether --verified_amplicon_fasta was supplied
    avg_verified_matches: Optional[float] = None
    avg_verified_variants: Optional[float] = None
    avg_verified_nonmatches: Optional[float] = None
    n_verified_matches_pairs: int = 0
    n_verified_variants_pairs: int = 0
    n_verified_nonmatches_pairs: int = 0
    verified_matches_by_position: list[PositionResult] = field(default_factory=list)
    verified_variants_by_position: list[PositionResult] = field(default_factory=list)
    verified_nonmatches_by_position: list[PositionResult] = field(default_factory=list)
    # "Variants" split further by whether the aligned monomer at that position is the
    # same substrate as the query's (a silent position, despite the polymer differing
    # elsewhere) or a genuinely different one (an actual substitution at this position).
    avg_matches_variants_substrate_match: Optional[float] = None
    avg_matches_variants_substrate_mismatch: Optional[float] = None
    n_matches_variants_substrate_match_pairs: int = 0
    n_matches_variants_substrate_mismatch_pairs: int = 0
    matches_variants_substrate_match_by_position: list[PositionResult] = field(default_factory=list)
    matches_variants_substrate_mismatch_by_position: list[PositionResult] = field(default_factory=list)
    avg_verified_variants_substrate_match: Optional[float] = None
    avg_verified_variants_substrate_mismatch: Optional[float] = None
    n_verified_variants_substrate_match_pairs: int = 0
    n_verified_variants_substrate_mismatch_pairs: int = 0
    verified_variants_substrate_match_by_position: list[PositionResult] = field(default_factory=list)
    verified_variants_substrate_mismatch_by_position: list[PositionResult] = field(default_factory=list)


def find_query_cluster(
    query_polymer: str,
    header_to_cluster: dict[str, int],
    header_to_polymer: dict[str, str],
) -> tuple[int, dict[int, list[str]]]:
    """Find the cluster id containing ``query_polymer`` and build a
    cluster_id -> [unique polymer strings] mapping for every cluster.

    :raises ValueError: if the query polymer isn't among the clustered unique
        polymers
    """
    cluster_to_polymers: dict[int, list[str]] = {}
    query_cluster_id: Optional[int] = None

    for header, polymer_string in header_to_polymer.items():
        cluster_id = header_to_cluster.get(header)
        if cluster_id is None:
            continue
        cluster_to_polymers.setdefault(cluster_id, []).append(polymer_string)
        if polymer_string == query_polymer:
            query_cluster_id = cluster_id

    if query_cluster_id is None:
        raise ValueError(
            f"Query polymer {query_polymer!r} was not found among the clustered "
            "unique polymers (check --separator and that the query matches the "
            "cluster FASTA's monomer naming exactly)"
        )

    return query_cluster_id, cluster_to_polymers


def occurrences_for_polymer(
    polymer_string: str,
    polymer_to_occurrences: dict[str, list[PredictedPolymer]],
) -> list[PredictedPolymer]:
    return polymer_to_occurrences.get(polymer_string, [])


# ---------------------------------------------------------------------------
# Non-match sampling
# ---------------------------------------------------------------------------

def sample_non_match_occurrence(
    cluster_id: int,
    cluster_to_polymers: dict[int, list[str]],
    polymer_to_occurrences: dict[str, list[PredictedPolymer]],
    amplicon_by_domain: dict[tuple[str, int, str, str], AmpliconRecord],
    rng: random.Random,
    max_attempts: int = 20,
) -> Optional[PredictedPolymer]:
    """Pick one random unique polymer from ``cluster_id`` and one random
    occurrence of it, to represent that cluster in the "non-match" baseline.

    Retries with different polymers/occurrences (up to ``max_attempts``) if a
    draw turns up an occurrence with no usable amplicon provenance at all
    (missing "Domains" data, or no amplicon actually resolves), since not
    every occurrence necessarily has one (e.g. missing PARAS/hmmscan
    results). Returns None if no usable occurrence was found after all
    attempts.
    """
    candidate_polymers = cluster_to_polymers.get(cluster_id, [])
    if not candidate_polymers:
        return None

    for _ in range(max_attempts):
        polymer_string = rng.choice(candidate_polymers)
        occurrences = occurrences_for_polymer(polymer_string, polymer_to_occurrences)
        if not occurrences:
            continue
        occurrence = rng.choice(occurrences)
        ordered = ordered_amplicons_for_occurrence(occurrence, amplicon_by_domain)
        if ordered is not None and any(amplicon is not None for amplicon in ordered):
            return occurrence

    return None


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

@dataclass
class _MatchGroup:
    """One unique "match" polymer (identical to, or containing, the query),
    and the offset at which the query starts within it -- e.g. offset=0 for
    an exact match, offset=2 if the query is a 13-mer starting at this
    member's 3rd monomer."""
    polymer: list[str]
    offset: int
    occurrences: list[PredictedPolymer] = field(default_factory=list)


@dataclass
class _VariantGroup:
    """One unique "variant" polymer, and its monomer-level alignment against
    the query (query_index -> variant_index pairs for every aligned, i.e.
    corresponding, column)."""
    polymer: list[str]
    position_pairs: list[tuple[int, int]]
    occurrences: list[PredictedPolymer] = field(default_factory=list)


def _summarize_by_position(
    query_monomers: list[str], by_position: dict[int, list[tuple[str, str]]]
) -> list[PositionResult]:
    results = [
        PositionResult(
            query_index=i,
            query_monomer=query_monomers[i],
            avg_identity=average_paired_identity(by_position.get(i, [])),
            n_pairs=len(by_position.get(i, [])),
        )
        for i in range(len(query_monomers))
    ]
    return results


def analyze_cluster(
    query_polymer: str,
    clusters_tsv: str,
    cluster_fasta: str,
    polymer_metadata_file: str,
    amplicon_fasta: str,
    separator: str = "|",
    seed: Optional[int] = 0,
    search_tsv: Optional[str] = None,
    verified_amplicon_fasta: Optional[str] = None,
) -> ClusterAnalysis:
    """Run the full analysis described in the module docstring.

    :param query_polymer: the query polymer, as a ``separator``-joined monomer
        string (must match one of the unique polymers that was clustered)
    :param clusters_tsv: path to NATU's clusters.tsv
    :param cluster_fasta: path to the FASTA of unique polymers that was
        clustered (header -> polymer sequence)
    :param polymer_metadata_file: path to the polymer metadata file with real
        genome occurrences (e.g. natu_all_polymers.txt), including a
        "Domains" column (see the module docstring)
    :param amplicon_fasta: path to the amplicon FASTA (post in-silico amplicon
        extraction)
    :param separator: monomer separator used throughout (default "|")
    :param seed: seed for the random representative draws used for the
        non-match baseline, for reproducibility. None for non-deterministic.
    :param search_tsv: optional path to a precomputed ``natu search`` output TSV
        (see ``load_search_alignments``) to reuse for variant/non-match
        monomer-level alignment instead of recomputing it. Any query/subject
        pair not found there falls back to the on-the-fly alignment as before.
    :param verified_amplicon_fasta: optional path to a FASTA of verified
        amplicons for the query's own cluster, one record per polymer
        position in polymer order (see ``load_verified_amplicons``). When
        given, adds a comparison of these against the in-silico amplicons of
        the match occurrences, variant occurrences, and non-match
        representatives (see the module docstring's step 7).
    """
    rng = random.Random(seed)

    search_alignments = load_search_alignments(search_tsv) if search_tsv else None

    header_to_cluster = load_clusters(clusters_tsv)
    header_to_polymer = parse_fasta_file(cluster_fasta)

    query_cluster_id, cluster_to_polymers = find_query_cluster(
        query_polymer, header_to_cluster, header_to_polymer
    )
    query_monomers = monomers_of(query_polymer, separator)

    verified_amplicons = (
        load_verified_amplicons(verified_amplicon_fasta, len(query_monomers))
        if verified_amplicon_fasta else None
    )

    result = ClusterAnalysis(
        query_polymer=query_polymer,
        query_monomers=query_monomers,
        cluster_id=query_cluster_id,
        cluster_members=cluster_to_polymers[query_cluster_id],
        used_search_tsv=search_tsv is not None,
        used_verified_amplicons=verified_amplicon_fasta is not None,
    )

    polymer_metadata = read_polymers(polymer_metadata_file, separator=separator)
    polymer_to_occurrences: dict[str, list[PredictedPolymer]] = {}
    for polymer in polymer_metadata:
        polymer_string = separator.join(polymer.polymer)
        polymer_to_occurrences.setdefault(polymer_string, []).append(polymer)

    amplicons = load_amplicons(amplicon_fasta, separator=separator)
    amplicon_by_domain = index_amplicons_by_domain(amplicons)

    def ordered_for(occurrence: PredictedPolymer) -> Optional[list[Optional[AmpliconRecord]]]:
        return ordered_amplicons_for_occurrence(occurrence, amplicon_by_domain)

    # Step 2 + 3: resolve cluster members to occurrences, split into match
    # groups (known offset, no alignment needed) and variant groups (aligned
    # once per unique polymer, reused for every occurrence of it).
    match_groups: list[_MatchGroup] = []
    variant_groups: list[_VariantGroup] = []

    for member_polymer in result.cluster_members:
        occurrences = occurrences_for_polymer(member_polymer, polymer_to_occurrences)
        if not occurrences:
            print(f"Warning: no genome occurrences found for cluster member {member_polymer!r}")
            continue

        member_monomers = monomers_of(member_polymer, separator)
        offset = match_offset(member_monomers, query_monomers)
        if offset is not None:
            match_groups.append(_MatchGroup(polymer=member_monomers, offset=offset, occurrences=occurrences))
            result.match_occurrences.extend(occurrences)
        else:
            position_pairs, from_search = resolve_variant_position_pairs(
                query_monomers, member_monomers, search_alignments
            )
            if from_search:
                result.n_position_pairs_from_search_tsv += 1
            else:
                result.n_position_pairs_computed += 1
            variant_groups.append(
                _VariantGroup(polymer=member_monomers, position_pairs=position_pairs, occurrences=occurrences)
            )
            result.variant_occurrences.extend(occurrences)

    # (occurrence, offset-of-query-within-it) for every match occurrence,
    # flattened across match groups -- the unit every comparison below is
    # built from.
    match_occurrences_with_offset: list[tuple[PredictedPolymer, int]] = [
        (occurrence, group.offset) for group in match_groups for occurrence in group.occurrences
    ]

    for occurrence, _ in match_occurrences_with_offset:
        if ordered_for(occurrence) is None:
            result.occurrences_without_provenance += 1
    for occurrence in result.variant_occurrences:
        if ordered_for(occurrence) is None:
            result.occurrences_without_provenance += 1

    # Step 4: position-paired amplicon similarity.

    # -- within matches: every pair of *distinct* match occurrences, compared
    # position-by-position at the query's own indices (offset already known
    # for each side).
    within_matches_by_position: dict[int, list[tuple[str, str]]] = {}
    for (occ_a, offset_a), (occ_b, offset_b) in combinations(match_occurrences_with_offset, 2):
        accumulate_by_query_position(
            within_matches_by_position,
            ordered_for(occ_a),
            ordered_for(occ_b),
            ((i, offset_a + i, offset_b + i) for i in range(len(query_monomers))),
        )

    # -- verified vs. matches: the verified reference against every match
    # occurrence's own in-silico amplicons, at the query's own indices.
    verified_matches_by_position: dict[int, list[tuple[str, str]]] = {}
    if verified_amplicons is not None:
        for match_occurrence, offset in match_occurrences_with_offset:
            accumulate_verified_by_position(
                verified_matches_by_position,
                verified_amplicons,
                ordered_for(match_occurrence),
                ((i, offset + i) for i in range(len(query_monomers))),
            )

    # -- matches vs. variants: each variant group's monomer-level alignment
    # (query_index -> variant_index) reused across all of that group's
    # occurrences, compared against every match occurrence. Also split into
    # substrate-match/substrate-mismatch subsets (see
    # ``split_position_pairs_by_substrate``) so a silent aligned position
    # doesn't get blended together with a genuine substitution.
    matches_variants_by_position: dict[int, list[tuple[str, str]]] = {}
    matches_variants_substrate_match_by_position: dict[int, list[tuple[str, str]]] = {}
    matches_variants_substrate_mismatch_by_position: dict[int, list[tuple[str, str]]] = {}
    verified_variants_by_position: dict[int, list[tuple[str, str]]] = {}
    verified_variants_substrate_match_by_position: dict[int, list[tuple[str, str]]] = {}
    verified_variants_substrate_mismatch_by_position: dict[int, list[tuple[str, str]]] = {}
    for group in variant_groups:
        if not group.position_pairs:
            continue
        substrate_match_pairs, substrate_mismatch_pairs = split_position_pairs_by_substrate(
            query_monomers, group.polymer, group.position_pairs
        )
        for variant_occurrence in group.occurrences:
            variant_ordered = ordered_for(variant_occurrence)
            if variant_ordered is None:
                continue
            if verified_amplicons is not None:
                accumulate_verified_by_position(
                    verified_variants_by_position, verified_amplicons, variant_ordered, group.position_pairs,
                )
                accumulate_verified_by_position(
                    verified_variants_substrate_match_by_position,
                    verified_amplicons, variant_ordered, substrate_match_pairs,
                )
                accumulate_verified_by_position(
                    verified_variants_substrate_mismatch_by_position,
                    verified_amplicons, variant_ordered, substrate_mismatch_pairs,
                )
            for match_occurrence, offset in match_occurrences_with_offset:
                match_ordered = ordered_for(match_occurrence)
                if match_ordered is None:
                    continue
                accumulate_by_query_position(
                    matches_variants_by_position,
                    match_ordered,
                    variant_ordered,
                    (
                        (query_index, offset + query_index, variant_index)
                        for query_index, variant_index in group.position_pairs
                    ),
                )
                accumulate_by_query_position(
                    matches_variants_substrate_match_by_position,
                    match_ordered,
                    variant_ordered,
                    (
                        (query_index, offset + query_index, variant_index)
                        for query_index, variant_index in substrate_match_pairs
                    ),
                )
                accumulate_by_query_position(
                    matches_variants_substrate_mismatch_by_position,
                    match_ordered,
                    variant_ordered,
                    (
                        (query_index, offset + query_index, variant_index)
                        for query_index, variant_index in substrate_mismatch_pairs
                    ),
                )

    # -- matches vs. non-matches: one representative occurrence per other
    # cluster, aligned to the query the same way a variant would be.
    matches_nonmatches_by_position: dict[int, list[tuple[str, str]]] = {}
    verified_nonmatches_by_position: dict[int, list[tuple[str, str]]] = {}
    other_cluster_ids = [cid for cid in cluster_to_polymers if cid != query_cluster_id]
    for cluster_id in other_cluster_ids:
        occurrence = sample_non_match_occurrence(
            cluster_id, cluster_to_polymers, polymer_to_occurrences, amplicon_by_domain, rng
        )
        if occurrence is None:
            result.skipped_non_match_clusters.append(cluster_id)
            continue
        result.non_match_cluster_ids.append(cluster_id)

        other_ordered = ordered_for(occurrence)
        position_pairs, from_search = resolve_variant_position_pairs(
            query_monomers, occurrence.polymer, search_alignments
        )
        if from_search:
            result.n_position_pairs_from_search_tsv += 1
        else:
            result.n_position_pairs_computed += 1
        if verified_amplicons is not None:
            accumulate_verified_by_position(
                verified_nonmatches_by_position, verified_amplicons, other_ordered, position_pairs,
            )
        for match_occurrence, offset in match_occurrences_with_offset:
            match_ordered = ordered_for(match_occurrence)
            if match_ordered is None:
                continue
            accumulate_by_query_position(
                matches_nonmatches_by_position,
                match_ordered,
                other_ordered,
                (
                    (query_index, offset + query_index, other_index)
                    for query_index, other_index in position_pairs
                ),
            )

    result.within_matches_by_position = _summarize_by_position(query_monomers, within_matches_by_position)
    result.matches_variants_by_position = _summarize_by_position(query_monomers, matches_variants_by_position)
    result.matches_variants_substrate_match_by_position = _summarize_by_position(
        query_monomers, matches_variants_substrate_match_by_position
    )
    result.matches_variants_substrate_mismatch_by_position = _summarize_by_position(
        query_monomers, matches_variants_substrate_mismatch_by_position
    )
    result.matches_nonmatches_by_position = _summarize_by_position(query_monomers, matches_nonmatches_by_position)
    result.verified_matches_by_position = _summarize_by_position(query_monomers, verified_matches_by_position)
    result.verified_variants_by_position = _summarize_by_position(query_monomers, verified_variants_by_position)
    result.verified_variants_substrate_match_by_position = _summarize_by_position(
        query_monomers, verified_variants_substrate_match_by_position
    )
    result.verified_variants_substrate_mismatch_by_position = _summarize_by_position(
        query_monomers, verified_variants_substrate_mismatch_by_position
    )
    result.verified_nonmatches_by_position = _summarize_by_position(query_monomers, verified_nonmatches_by_position)

    all_within_matches_pairs = [pair for pairs in within_matches_by_position.values() for pair in pairs]
    all_matches_variants_pairs = [pair for pairs in matches_variants_by_position.values() for pair in pairs]
    all_matches_variants_substrate_match_pairs = [
        pair for pairs in matches_variants_substrate_match_by_position.values() for pair in pairs
    ]
    all_matches_variants_substrate_mismatch_pairs = [
        pair for pairs in matches_variants_substrate_mismatch_by_position.values() for pair in pairs
    ]
    all_matches_nonmatches_pairs = [pair for pairs in matches_nonmatches_by_position.values() for pair in pairs]
    all_verified_matches_pairs = [pair for pairs in verified_matches_by_position.values() for pair in pairs]
    all_verified_variants_pairs = [pair for pairs in verified_variants_by_position.values() for pair in pairs]
    all_verified_variants_substrate_match_pairs = [
        pair for pairs in verified_variants_substrate_match_by_position.values() for pair in pairs
    ]
    all_verified_variants_substrate_mismatch_pairs = [
        pair for pairs in verified_variants_substrate_mismatch_by_position.values() for pair in pairs
    ]
    all_verified_nonmatches_pairs = [pair for pairs in verified_nonmatches_by_position.values() for pair in pairs]

    result.avg_within_matches = average_paired_identity(all_within_matches_pairs)
    result.n_within_matches_pairs = len(all_within_matches_pairs)

    result.avg_matches_variants = average_paired_identity(all_matches_variants_pairs)
    result.n_matches_variants_pairs = len(all_matches_variants_pairs)

    result.avg_matches_variants_substrate_match = average_paired_identity(all_matches_variants_substrate_match_pairs)
    result.n_matches_variants_substrate_match_pairs = len(all_matches_variants_substrate_match_pairs)

    result.avg_matches_variants_substrate_mismatch = average_paired_identity(all_matches_variants_substrate_mismatch_pairs)
    result.n_matches_variants_substrate_mismatch_pairs = len(all_matches_variants_substrate_mismatch_pairs)

    result.avg_matches_nonmatches = average_paired_identity(all_matches_nonmatches_pairs)
    result.n_matches_nonmatches_pairs = len(all_matches_nonmatches_pairs)

    result.avg_verified_matches = average_paired_identity(all_verified_matches_pairs)
    result.n_verified_matches_pairs = len(all_verified_matches_pairs)

    result.avg_verified_variants = average_paired_identity(all_verified_variants_pairs)
    result.n_verified_variants_pairs = len(all_verified_variants_pairs)

    result.avg_verified_variants_substrate_match = average_paired_identity(all_verified_variants_substrate_match_pairs)
    result.n_verified_variants_substrate_match_pairs = len(all_verified_variants_substrate_match_pairs)

    result.avg_verified_variants_substrate_mismatch = average_paired_identity(all_verified_variants_substrate_mismatch_pairs)
    result.n_verified_variants_substrate_mismatch_pairs = len(all_verified_variants_substrate_mismatch_pairs)

    result.avg_verified_nonmatches = average_paired_identity(all_verified_nonmatches_pairs)
    result.n_verified_nonmatches_pairs = len(all_verified_nonmatches_pairs)

    return result


def suggest_variant_range(result: ClusterAnalysis) -> Optional[tuple[float, float]]:
    """Heuristic suggested percent-identity range for calling something a
    "predicted variant" of the query polymer, based on where the observed
    matches-vs-variants similarity sits between the matches-vs-nonmatches
    ("unrelated") baseline and the within-matches ("same") baseline.

    The lower bound is the midpoint between the non-match baseline and the
    observed matches-variants average; the upper bound is the midpoint between
    the matches-variants average and the within-matches baseline. This is a
    single-cluster heuristic, not a validated cutoff -- check it against a few
    more clusters before relying on it.
    """
    if None in (result.avg_within_matches, result.avg_matches_variants, result.avg_matches_nonmatches):
        return None

    lower = (result.avg_matches_nonmatches + result.avg_matches_variants) / 2
    upper = (result.avg_matches_variants + result.avg_within_matches) / 2
    return lower, upper


def suggest_variant_range_verified(result: ClusterAnalysis) -> Optional[tuple[float, float]]:
    """Like ``suggest_variant_range``, but anchored on the verified reference
    (see ``--verified_amplicon_fasta``) instead of the within-matches
    cross-occurrence comparison -- i.e. using how well the in-silico
    amplicons of matches/variants/non-matches agree with a trusted reference,
    rather than how similar those groups' in-silico amplicons are to each
    other. Returns None if ``--verified_amplicon_fasta`` wasn't supplied, or
    if any of the three baselines had no comparable pairs.
    """
    if None in (result.avg_verified_matches, result.avg_verified_variants, result.avg_verified_nonmatches):
        return None

    lower = (result.avg_verified_nonmatches + result.avg_verified_variants) / 2
    upper = (result.avg_verified_variants + result.avg_verified_matches) / 2
    return lower, upper


def print_report(result: ClusterAnalysis) -> None:
    print(f"Query polymer: {result.query_polymer}")
    print(f"Cluster id: {result.cluster_id} ({len(result.cluster_members)} unique polymer members)")
    print(f"Match occurrences: {len(result.match_occurrences)}")
    print(f"Variant occurrences: {len(result.variant_occurrences)}")
    if result.occurrences_without_provenance:
        print(
            f"Warning: {result.occurrences_without_provenance} match/variant occurrence(s) had no "
            "domain provenance (\"Domains\" column) and were skipped from every comparison below -- "
            "re-run get_polymers on a version that writes it if this number looks large."
        )
    if result.used_search_tsv:
        total_alignments = result.n_position_pairs_from_search_tsv + result.n_position_pairs_computed
        print(
            f"Monomer alignments (variants + non-match representatives): "
            f"{result.n_position_pairs_from_search_tsv}/{total_alignments} reused from --search_tsv, "
            f"{result.n_position_pairs_computed} computed on the fly (pair not found there)"
        )
    print()

    def fmt(avg: Optional[float], n_pairs: int) -> str:
        return f"{avg:.2f}% (n={n_pairs} position-paired comparisons)" if avg is not None else f"n/a (n={n_pairs})"

    print(f"Avg amplicon identity within matches:            {fmt(result.avg_within_matches, result.n_within_matches_pairs)}")
    print(f"Avg amplicon identity matches vs. variants:      {fmt(result.avg_matches_variants, result.n_matches_variants_pairs)}")
    print(f"  -- substrate match (same monomer at that position):          "
          f"{fmt(result.avg_matches_variants_substrate_match, result.n_matches_variants_substrate_match_pairs)}")
    print(f"  -- substrate mismatch (different monomer at that position):  "
          f"{fmt(result.avg_matches_variants_substrate_mismatch, result.n_matches_variants_substrate_mismatch_pairs)}")
    print(f"Avg amplicon identity matches vs. non-matches:   {fmt(result.avg_matches_nonmatches, result.n_matches_nonmatches_pairs)}")
    print(f"  (non-match baseline from {len(result.non_match_cluster_ids)} other cluster(s); "
          f"{len(result.skipped_non_match_clusters)} cluster(s) skipped, no usable amplicon found)")
    print()

    variant_range = suggest_variant_range(result)
    if variant_range is not None:
        print(f"Suggested predicted-variant identity range: {variant_range[0]:.2f}% - {variant_range[1]:.2f}%")
        print("(heuristic based on this cluster only -- validate against more clusters before adopting as a general cutoff)")
    else:
        print("Could not suggest a variant range: one or more baselines had no comparable amplicon pairs.")
    print()

    print("Per-position amplicon identity (query monomer -> within matches / "
          "vs. variants [substrate match / substrate mismatch] / vs. non-matches):")
    for within, variants_match, variants_mismatch, nonmatches in zip(
        result.within_matches_by_position,
        result.matches_variants_substrate_match_by_position,
        result.matches_variants_substrate_mismatch_by_position,
        result.matches_nonmatches_by_position,
    ):

        def cell(position_result: PositionResult) -> str:
            if position_result.avg_identity is None:
                return f"n/a (n={position_result.n_pairs})"
            return f"{position_result.avg_identity:5.1f}% (n={position_result.n_pairs})"

        print(
            f"  [{within.query_index:2d}] {within.query_monomer:<28s} "
            f"within={cell(within)}  vs_variants[match]={cell(variants_match)}  "
            f"vs_variants[mismatch]={cell(variants_mismatch)}  vs_nonmatches={cell(nonmatches)}"
        )

    if result.used_verified_amplicons:
        print()
        print(f"Avg amplicon identity, verified vs. matches:      {fmt(result.avg_verified_matches, result.n_verified_matches_pairs)}")
        print(f"Avg amplicon identity, verified vs. variants:     {fmt(result.avg_verified_variants, result.n_verified_variants_pairs)}")
        print(f"  -- substrate match (same monomer at that position):          "
              f"{fmt(result.avg_verified_variants_substrate_match, result.n_verified_variants_substrate_match_pairs)}")
        print(f"  -- substrate mismatch (different monomer at that position):  "
              f"{fmt(result.avg_verified_variants_substrate_mismatch, result.n_verified_variants_substrate_mismatch_pairs)}")
        print(f"Avg amplicon identity, verified vs. non-matches:  {fmt(result.avg_verified_nonmatches, result.n_verified_nonmatches_pairs)}")
        print()

        verified_range = suggest_variant_range_verified(result)
        if verified_range is not None:
            print(f"Suggested predicted-variant identity range (verified-anchored): {verified_range[0]:.2f}% - {verified_range[1]:.2f}%")
            print("(anchored on --verified_amplicon_fasta rather than cross-occurrence comparison)")
        else:
            print("Could not suggest a verified-anchored variant range: one or more baselines had no comparable amplicon pairs.")
        print()

        print("Per-position amplicon identity vs. verified reference "
              "(query monomer -> vs. matches / vs. variants [substrate match / substrate mismatch] / vs. non-matches):")
        for verified_m, verified_v_match, verified_v_mismatch, verified_n in zip(
            result.verified_matches_by_position,
            result.verified_variants_substrate_match_by_position,
            result.verified_variants_substrate_mismatch_by_position,
            result.verified_nonmatches_by_position,
        ):

            def cell(position_result: PositionResult) -> str:
                if position_result.avg_identity is None:
                    return f"n/a (n={position_result.n_pairs})"
                return f"{position_result.avg_identity:5.1f}% (n={position_result.n_pairs})"

            print(
                f"  [{verified_m.query_index:2d}] {verified_m.query_monomer:<28s} "
                f"vs_matches={cell(verified_m)}  vs_variants[match]={cell(verified_v_match)}  "
                f"vs_variants[mismatch]={cell(verified_v_mismatch)}  vs_nonmatches={cell(verified_n)}"
            )


def write_occurrence_table(occurrences: list[PredictedPolymer], label: str, out_file: str, separator: str = "|") -> None:
    with open(out_file, "w") as out:
        out.write("Accession\tRegion\tCandidate\tLocus tag\tPolymer\tGroup\n")
        for occurrence in occurrences:
            out.write(
                f"{occurrence.accession}\t{occurrence.region}\t"
                f"{occurrence.candidate_cluster or ''}\t{occurrence.locus_tag or ''}\t"
                f"{separator.join(occurrence.polymer)}\t{label}\n"
            )


def write_ranges_table(result: ClusterAnalysis, out_file: str) -> None:
    """Write the aggregate identity averages and the suggested predicted-variant
    identity range(s) to a TSV (columns: metric, avg_identity, n_pairs), rather
    than leaving them only in the printed report -- these are usually the
    numbers someone wants to gather across many clusters, or hand to another
    script, rather than read off the console one cluster at a time.

    ``avg_identity``/``n_pairs`` are "NA" for a metric that had no comparable
    pairs. The two suggested-range rows have "NA" for ``n_pairs`` -- a range
    is an interval derived from the metrics above it, not itself an average
    over a set of pairs. The verified-vs-* rows and the verified-anchored
    range are only written when ``--verified_amplicon_fasta`` was supplied
    (``result.used_verified_amplicons``).
    """
    def fmt_avg(avg: Optional[float]) -> str:
        return f"{avg:.4f}" if avg is not None else "NA"

    rows: list[tuple[str, str, str]] = [
        ("within_matches", fmt_avg(result.avg_within_matches), str(result.n_within_matches_pairs)),
        ("matches_variants", fmt_avg(result.avg_matches_variants), str(result.n_matches_variants_pairs)),
        (
            "matches_variants_substrate_match",
            fmt_avg(result.avg_matches_variants_substrate_match),
            str(result.n_matches_variants_substrate_match_pairs),
        ),
        (
            "matches_variants_substrate_mismatch",
            fmt_avg(result.avg_matches_variants_substrate_mismatch),
            str(result.n_matches_variants_substrate_mismatch_pairs),
        ),
        ("matches_nonmatches", fmt_avg(result.avg_matches_nonmatches), str(result.n_matches_nonmatches_pairs)),
    ]

    variant_range = suggest_variant_range(result)
    rows.append(("suggested_variant_range_lower", f"{variant_range[0]:.4f}" if variant_range else "NA", "NA"))
    rows.append(("suggested_variant_range_upper", f"{variant_range[1]:.4f}" if variant_range else "NA", "NA"))

    if result.used_verified_amplicons:
        rows.append(("verified_matches", fmt_avg(result.avg_verified_matches), str(result.n_verified_matches_pairs)))
        rows.append(("verified_variants", fmt_avg(result.avg_verified_variants), str(result.n_verified_variants_pairs)))
        rows.append(
            (
                "verified_variants_substrate_match",
                fmt_avg(result.avg_verified_variants_substrate_match),
                str(result.n_verified_variants_substrate_match_pairs),
            )
        )
        rows.append(
            (
                "verified_variants_substrate_mismatch",
                fmt_avg(result.avg_verified_variants_substrate_mismatch),
                str(result.n_verified_variants_substrate_mismatch_pairs),
            )
        )
        rows.append(
            ("verified_nonmatches", fmt_avg(result.avg_verified_nonmatches), str(result.n_verified_nonmatches_pairs))
        )

        verified_range = suggest_variant_range_verified(result)
        rows.append(
            ("suggested_variant_range_verified_lower", f"{verified_range[0]:.4f}" if verified_range else "NA", "NA")
        )
        rows.append(
            ("suggested_variant_range_verified_upper", f"{verified_range[1]:.4f}" if verified_range else "NA", "NA")
        )

    with open(out_file, "w") as out:
        out.write("metric\tavg_identity\tn_pairs\n")
        for metric, avg, n_pairs in rows:
            out.write(f"{metric}\t{avg}\t{n_pairs}\n")


def write_per_position_table(result: ClusterAnalysis, out_file: str) -> None:
    """Write the per-position identity breakdown -- the data the ranges in
    ``write_ranges_table`` are derived from -- to a TSV, one row per query
    polymer position, for archiving or plotting across many clusters.

    The ``verified_*`` columns are always present, for a consistent schema
    across runs regardless of whether any given run used
    ``--verified_amplicon_fasta``; they're filled with "NA" when it wasn't.
    """
    def fmt(position_result: PositionResult) -> tuple[str, str]:
        avg = f"{position_result.avg_identity:.4f}" if position_result.avg_identity is not None else "NA"
        return avg, str(position_result.n_pairs)

    header = [
        "query_index", "query_monomer",
        "within_matches_avg", "within_matches_n",
        "matches_variants_substrate_match_avg", "matches_variants_substrate_match_n",
        "matches_variants_substrate_mismatch_avg", "matches_variants_substrate_mismatch_n",
        "matches_nonmatches_avg", "matches_nonmatches_n",
        "verified_matches_avg", "verified_matches_n",
        "verified_variants_substrate_match_avg", "verified_variants_substrate_match_n",
        "verified_variants_substrate_mismatch_avg", "verified_variants_substrate_mismatch_n",
        "verified_nonmatches_avg", "verified_nonmatches_n",
    ]

    with open(out_file, "w") as out:
        out.write("\t".join(header) + "\n")
        for i in range(len(result.query_monomers)):
            within_avg, within_n = fmt(result.within_matches_by_position[i])
            vm_avg, vm_n = fmt(result.matches_variants_substrate_match_by_position[i])
            vmm_avg, vmm_n = fmt(result.matches_variants_substrate_mismatch_by_position[i])
            nm_avg, nm_n = fmt(result.matches_nonmatches_by_position[i])
            if result.used_verified_amplicons:
                verm_avg, verm_n = fmt(result.verified_matches_by_position[i])
                vervm_avg, vervm_n = fmt(result.verified_variants_substrate_match_by_position[i])
                vervmm_avg, vervmm_n = fmt(result.verified_variants_substrate_mismatch_by_position[i])
                vernm_avg, vernm_n = fmt(result.verified_nonmatches_by_position[i])
            else:
                verm_avg = verm_n = vervm_avg = vervm_n = vervmm_avg = vervmm_n = vernm_avg = vernm_n = "NA"
            row = [
                str(i), result.query_monomers[i],
                within_avg, within_n,
                vm_avg, vm_n,
                vmm_avg, vmm_n,
                nm_avg, nm_n,
                verm_avg, verm_n,
                vervm_avg, vervm_n,
                vervmm_avg, vervmm_n,
                vernm_avg, vernm_n,
            ]
            out.write("\t".join(row) + "\n")


def parse_arguments() -> Namespace:
    parser = ArgumentParser(
        description="Link a NATU cluster to real genome occurrences, split them into "
        "matches/variants of a query polymer, and compare AMP-binding domain amplicon "
        "similarity within/between those groups and other clusters, position by position."
    )
    parser.add_argument("-q", "--query_polymer", required=True, type=str,
                        help="Query polymer as a separator-joined monomer string")
    parser.add_argument("-c", "--clusters", required=True, type=str,
                        help="Path to NATU's clusters.tsv")
    parser.add_argument("-f", "--cluster_fasta", required=True, type=str,
                        help="Path to the FASTA file of unique polymers that was clustered")
    parser.add_argument("-m", "--polymer_metadata", required=True, type=str,
                        help="Path to the polymer metadata file with real genome occurrences "
                        "(e.g. natu_all_polymers.txt), including a Domains column")
    parser.add_argument("-a", "--amplicon_fasta", required=True, type=str,
                        help="Path to the amplicon FASTA (after in-silico amplicon extraction)")
    parser.add_argument("--search_tsv", type=str, default=None,
                        help="Optional path to a precomputed 'natu search' output TSV "
                        "(columns: query, subject, bitscore, aligned_q, s_aligned_s) to reuse "
                        "for variant/non-match monomer-level alignment instead of recomputing "
                        "it here. For guaranteed coverage of every cluster member and any "
                        "randomly-sampled non-match representative, generate it with "
                        "'natu search --threshold 0' against the same --cluster_fasta, ideally "
                        "with the same substitution matrix/config used for the 'natu cluster' "
                        "run that produced --clusters. Any pair not found in this file falls "
                        "back to this script's own on-the-fly alignment.")
    parser.add_argument("--verified_amplicon_fasta", type=str, default=None,
                        help="Optional path to a FASTA of verified amplicons for the query's "
                        "own cluster: one record per polymer position, in polymer order "
                        "(position 0 first; headers aren't parsed). When given, the report "
                        "adds a comparison of these against the in-silico amplicons of the "
                        "match occurrences, variant occurrences, and non-match representatives "
                        "-- i.e. how well the automated extraction agrees with a trusted "
                        "reference, rather than how similar those groups are to each other. "
                        "Record count must equal the query polymer's monomer count.")
    parser.add_argument("-s", "--separator", default="|", type=str,
                        help="Monomer separator (default: '|')")
    parser.add_argument("--seed", default=0, type=int,
                        help="Random seed for non-match sampling (default: 0). Use a negative "
                        "value for a non-deterministic seed.")
    parser.add_argument("-o", "--out_dir", type=str, default=None,
                        help="If given, write match/variant occurrence tables, the aggregate "
                        "identity averages and suggested variant range(s) (ranges.tsv), and the "
                        "per-position identity breakdown (per_position.tsv) to this directory")

    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    seed = None if args.seed < 0 else args.seed

    result = analyze_cluster(
        query_polymer=args.query_polymer,
        clusters_tsv=args.clusters,
        cluster_fasta=args.cluster_fasta,
        polymer_metadata_file=args.polymer_metadata,
        amplicon_fasta=args.amplicon_fasta,
        separator=args.separator,
        seed=seed,
        search_tsv=args.search_tsv,
        verified_amplicon_fasta=args.verified_amplicon_fasta,
    )

    print_report(result)

    if args.out_dir:
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        write_occurrence_table(result.match_occurrences, "match", str(out_dir / "matches.tsv"), args.separator)
        write_occurrence_table(result.variant_occurrences, "variant", str(out_dir / "variants.tsv"), args.separator)
        write_ranges_table(result, str(out_dir / "ranges.tsv"))
        write_per_position_table(result, str(out_dir / "per_position.tsv"))
        print(f"\nWrote matches.tsv, variants.tsv, ranges.tsv, and per_position.tsv to {out_dir}")


if __name__ == "__main__":
    main()
