"""
Extract antiSMASH AMP-binding (adenylation) domain sequences from a folder of
antiSMASH JSON files and write them all to a single FASTA file.

Each sequence is named as a single whitespace-free token:

    <accession>|region_<region>|candidate_<candidates>|<locus_tag>|<domain_index>

which mirrors the Accession / Region / Candidate columns used by
`poliwag.antismash_processing.get_polymers.write_polymers` (e.g. the
`all_polymers.txt` / `unique_polymers.txt` files that script produces), so a
domain here can be matched to its polymer row on Accession + Region +
Candidate. `<candidates>` is `-`-joined when a domain falls in more than one
overlapping (interleaved) candidate cluster, and `NA` if it falls in none.

Unlike `poliwag.antismash_processing.get_a_domains`, this script expects a
flat folder of .json files (one antiSMASH run per file), not a folder of
per-genome subfolders.
"""

import json
import traceback
from argparse import ArgumentParser, Namespace
from dataclasses import dataclass
from typing import Any, Optional, TextIO

from poliwag.antismash_processing.antismash_json import (
    domain_suffix,
    get_cand_clusters_from_location,
    get_cds_from_name,
    parse_location,
    region_from_cds,
)
from poliwag.antismash_processing.location import Location
from poliwag.utils import iterate_over_dir
from tqdm import tqdm


@dataclass
class AmpBindingDomain:
    """A single AMP-binding (adenylation) domain extracted from antiSMASH JSON."""

    accession: str
    region: int
    candidate_clusters: list[int]
    locus_tag: str
    domain_index: str
    location: Location
    protein_start: int
    protein_end: int
    sequence: str
    paras_substrate: Optional[str] = None
    paras_confidence: Optional[float] = None
    evalue: Optional[str] = None
    score: Optional[str] = None

    @property
    def candidate_label(self) -> str:
        if self.candidate_clusters:
            return "-".join(str(c) for c in sorted(self.candidate_clusters))
        return "NA"

    @property
    def fasta_id(self) -> str:
        """Single whitespace-free token that cross-links to a polymer table
        (see module docstring)."""
        return (
            f"{self.accession}|region_{self.region}|"
            f"candidate_{self.candidate_label}|{self.locus_tag}|{self.domain_index}"
        )

    def write_to_fasta(self, fasta_file: TextIO) -> None:
        confidence = f"{self.paras_confidence:.3f}" if self.paras_confidence is not None else "NA"
        description = (
            f"locus_tag={self.locus_tag}"
            f"\tnt_coordinates={self.location}"
            f"\taa_coordinates={self.protein_start}-{self.protein_end}"
            f"\tparas_prediction={self.paras_substrate or 'NA'}"
            f"\tparas_confidence={confidence}"
        )
        if self.evalue is not None:
            description += f"\tevalue={self.evalue}"
        if self.score is not None:
            description += f"\tscore={self.score}"
        fasta_file.write(f">{self.fasta_id}\t{description}\n{self.sequence}\n")


def _get_paras_prediction(domain_id: str, record: dict[str, Any]) -> tuple[Optional[str], Optional[float]]:
    """Look up the PARAS substrate prediction for a domain from the record's
    nrps_pks module results.

    :param domain_id: the aSDomain's domain_id, e.g. '..._AMP-binding.1'
    :param record: antiSMASH record in JSON format
    :return: (predicted_substrate, confidence), or (None, None) if PARAS wasn't run
        for this domain (e.g. no NRPS/PKS module predictions in this record)
    """
    domain_predictions = (
        record.get("modules", {})
        .get("antismash.modules.nrps_pks", {})
        .get("domain_predictions", {})
    )
    paras_result = domain_predictions.get(domain_id, {}).get("paras")
    if not paras_result:
        return None, None
    return paras_result.get("predicted_substrate"), paras_result.get("confidence")


def amp_binding_domains_from_json(antismash_json_file: str) -> list[AmpBindingDomain]:
    """Extract all AMP-binding domains from a single antiSMASH JSON file.

    :param antismash_json_file: path to an antiSMASH JSON output file
    :return: list of AMP-binding domains found in the file
    """
    with open(antismash_json_file) as f:
        antismash_data = json.load(f)

    domains: list[AmpBindingDomain] = []
    for record in antismash_data["records"]:
        accession = record["id"]
        for feature in record["features"]:
            if feature["type"] != "aSDomain":
                continue

            qualifiers: dict[str, Any] = feature["qualifiers"]
            domain_ids = qualifiers.get("domain_id")
            if not domain_ids or "AMP-binding" not in domain_ids[0]:
                continue

            domain_id = domain_ids[0]
            locus_tag = qualifiers["locus_tag"][0]
            location = parse_location(feature["location"])

            cds = get_cds_from_name(locus_tag, record)
            if cds is None:
                print(f"Warning: could not find CDS for {locus_tag} in {accession}, skipping")
                continue

            region = region_from_cds(cds, record)
            candidate_clusters = get_cand_clusters_from_location(location, record)
            paras_substrate, paras_confidence = _get_paras_prediction(domain_id, record)

            domains.append(
                AmpBindingDomain(
                    accession=accession,
                    region=region,
                    candidate_clusters=candidate_clusters,
                    locus_tag=locus_tag,
                    domain_index=domain_suffix(domain_id),
                    location=location,
                    protein_start=int(qualifiers["protein_start"][0]),
                    protein_end=int(qualifiers["protein_end"][0]),
                    sequence=qualifiers["translation"][0],
                    paras_substrate=paras_substrate,
                    paras_confidence=paras_confidence,
                    evalue=qualifiers.get("evalue", [None])[0],
                    score=qualifiers.get("score", [None])[0],
                )
            )
    return domains


def amp_binding_domains_from_folder(json_folder: str) -> list[AmpBindingDomain]:
    """Extract AMP-binding domains from every .json file in a flat folder.

    :param json_folder: folder containing antiSMASH .json output files
    :return: list of AMP-binding domains found across all files
    """
    all_domains: list[AmpBindingDomain] = []
    json_files = list(iterate_over_dir(json_folder, ".json"))
    for file_name, file_path in tqdm(json_files, desc="Extracting AMP-binding domains", unit="file"):
        try:
            domains = amp_binding_domains_from_json(file_path)
            all_domains.extend(domains)
        except Exception:
            print(f"Could not process {file_name}")
            print(traceback.format_exc())
    return all_domains


def write_fasta(domains: list[AmpBindingDomain], out_file: str) -> None:
    """Write a list of AMP-binding domains to a FASTA file."""
    with open(out_file, "w") as out:
        for domain in domains:
            domain.write_to_fasta(out)


def parse_arguments() -> Namespace:
    parser = ArgumentParser(
        description="Extract AMP-binding domain sequences from a folder of antiSMASH "
        "JSON files into a single FASTA file, named for cross-linking to a polymer file."
    )
    parser.add_argument(
        "-i", "--json_folder", type=str, required=True,
        help="Folder containing antiSMASH .json output files (one file per genome)",
    )
    parser.add_argument(
        "-o", "--out_file", type=str, required=True,
        help="Path to output FASTA file",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    domains = amp_binding_domains_from_folder(args.json_folder)
    write_fasta(domains, args.out_file)
    print(f"Wrote {len(domains)} AMP-binding domain sequences to {args.out_file}")


if __name__ == "__main__":
    main()
