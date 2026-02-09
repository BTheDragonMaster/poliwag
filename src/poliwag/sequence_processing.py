from typing import Optional

from Bio.SeqRecord import SeqRecord


def get_organism(record: SeqRecord) -> str:
    organism_name = record.annotations.get("organism")
    return organism_name


def get_assembly_accession(record: SeqRecord) -> Optional[str]:
    dbxrefs = record.annotations.get("dbxrefs", [])

    assembly_accession = None
    for entry in dbxrefs:
        if entry.startswith("Assembly:"):
            assembly_accession = entry.split("Assembly:")[1].strip()
            break

    return assembly_accession
