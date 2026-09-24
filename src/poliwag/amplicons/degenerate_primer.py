"""
Shared helpers for turning a degenerate DNA primer into amino-acid-level
information: which amino acids each of its (possibly ambiguous) codons can
encode, and from that either an exact amino-acid regex or a single "most
common" consensus peptide.

Used by both ``degenerate_primer_amplicons.py`` (exact regex matching) and
``in_silico_amplicon_sequencing.py`` (local pairwise alignment to a
consensus peptide), so the codon-expansion logic lives in one place.
"""

from collections import Counter
from dataclasses import dataclass
from itertools import product
from typing import List, Tuple

from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Seq import Seq


def split_into_codons(primer: str) -> Tuple[List[str], str]:
    """Split a primer into its complete codons, plus any trailing partial codon.

    :param primer: DNA primer sequence (IUPAC ambiguity codes allowed)
    :return: (list of 3-base codons, leftover bases that don't complete a codon)
    """
    n_complete = len(primer) - (len(primer) % 3)
    codons = [primer[i:i + 3] for i in range(0, n_complete, 3)]
    return codons, primer[n_complete:]


def expand_iupac_codon(codon: str) -> List[str]:
    """Expand a (possibly degenerate) codon into every concrete codon it represents.

    :param codon: 3-base codon, IUPAC ambiguity codes allowed
    :return: every concrete A/C/G/T codon consistent with ``codon``
    """
    choices = [ambiguous_dna_values[base.upper()] for base in codon]
    return ["".join(bases) for bases in product(*choices)]


def codon_amino_acid_counts(codon: str) -> Counter:
    """Translate every concrete codon a degenerate codon represents, and
    count how many of them encode each amino acid.

    :param codon: 3-base codon, IUPAC ambiguity codes allowed
    :return: Counter mapping one-letter amino acid (``*`` for stop) to the
        number of concrete codons that translate to it
    """
    return Counter(str(Seq(concrete).translate()) for concrete in expand_iupac_codon(codon))


@dataclass
class CodonTranslation:
    """The amino acid(s) a single (possibly degenerate) primer codon can encode."""

    codon: str
    amino_acids: Counter  # amino acid -> number of concrete codons encoding it

    @property
    def possible(self) -> List[str]:
        """Every amino acid this codon could encode, most common first."""
        return [aa for aa, _count in self.amino_acids.most_common()]

    @property
    def is_tie(self) -> bool:
        """Whether more than one amino acid is (jointly) the most common."""
        best_count = self.amino_acids.most_common(1)[0][1]
        return sum(1 for count in self.amino_acids.values() if count == best_count) > 1

    @property
    def consensus(self) -> str:
        """The single most common amino acid this codon encodes.

        Ties (more than one amino acid encoded by the same number of
        concrete codons) are broken alphabetically, so this is one of
        possibly several equally "most common" amino acids -- a
        deterministic but otherwise arbitrary choice among them.
        """
        best_count = self.amino_acids.most_common(1)[0][1]
        tied = sorted(aa for aa, count in self.amino_acids.items() if count == best_count)
        return tied[0]


def translate_codons(primer: str) -> Tuple[List[CodonTranslation], str]:
    """Translate every complete codon in a (possibly degenerate) primer.

    :param primer: DNA primer sequence (IUPAC ambiguity codes allowed)
    :return: (one ``CodonTranslation`` per complete codon, leftover bases
        that don't complete a codon and so aren't translated)
    """
    codons, leftover = split_into_codons(primer)
    translations = [CodonTranslation(codon, codon_amino_acid_counts(codon)) for codon in codons]
    return translations, leftover
