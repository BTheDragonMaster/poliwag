from dataclasses import dataclass
from typing import IO

from poliwag.antismash_processing.location import Location


@dataclass
class AdenylationDomain:
    sequence: str
    location: Location
    region: int
    candidate_clusters: list[int]
    locus_tag: str
    protein_start: int
    protein_end: int
    paras_prediction: str | None
    paras_confidence: float | None

    def write_to_fasta(self, fasta_file: IO):
        fasta_file.write(f">{self.locus_tag}\t{str(self.location)}\tprot({self.protein_start}-{self.protein_end})\tregion_{self.region}\tcand_clusters_{','.join(list(map(str, self.candidate_clusters)))}\t{self.paras_prediction}\t{self.paras_confidence}\n{self.sequence}\n")