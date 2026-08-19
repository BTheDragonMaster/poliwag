"""
Class for storing a predicted polymer
"""

from typing import Optional


class PredictedPolymer:
    def __init__(self, polymer: list[str],
                 accession: str,
                 region_number: int,
                 candidate_cluster: Optional[int] = None,
                 locus_tag: Optional[str] = None):
        self.polymer = polymer
        self.accession = accession
        self.region = region_number
        self.candidate_cluster = candidate_cluster
        self.locus_tag = locus_tag

    def __repr__(self):
        return ' | '.join(self.polymer)

    def __str__(self):
        return ' | '.join(self.polymer)