"""
Class for storing a predicted polymer
"""

from typing import Optional


class PredictedPolymer:
    def __init__(self, polymer: list[str],
                 accession: str,
                 region_number: int,
                 candidate_cluster: Optional[int] = None,
                 locus_tag: Optional[str] = None,
                 domain_ids: Optional[list[Optional[tuple[str, str]]]] = None):
        self.polymer = polymer
        self.accession = accession
        self.region = region_number
        self.candidate_cluster = candidate_cluster
        self.locus_tag = locus_tag
        # Parallel to `polymer`: for each monomer, the (locus_tag, domain_index)
        # of the AMP-binding domain that produced it (matching the identifiers
        # used in amp_binding_domains_to_fasta.py's AmpBindingDomain.fasta_id),
        # or None for a monomer with no corresponding AMP-binding domain (e.g.
        # a PKS_AT/CAL-derived building block). None (rather than a list) means
        # this polymer's provenance wasn't recorded at all (e.g. an older file
        # read without a "Domains" column).
        self.domain_ids = domain_ids

    def __repr__(self):
        return '|'.join(self.polymer)

    def __str__(self):
        return '|'.join(self.polymer)