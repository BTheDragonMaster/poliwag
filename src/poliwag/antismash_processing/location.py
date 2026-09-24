from dataclasses import dataclass


@dataclass
class Location:
    """
    Class to store a genomic location
    """
    start: int
    end: int
    strand: str

    def is_contained_in(self, location: "Location", strand_must_match: bool = False) -> bool:
        if strand_must_match and self.strand != location.strand:
            return False
        else:
            if self.start >= location.start and self.end <= location.end:
                return True
            else:
                return False

    def __str__(self):
        return f"{self.start}-{self.end}({self.strand})"