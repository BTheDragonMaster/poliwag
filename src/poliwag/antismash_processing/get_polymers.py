"""
Get predicted NRPS and PKS polymers from antiSMASH output
"""

import os
from argparse import ArgumentParser, Namespace
import traceback

from poliwag.utils import iterate_over_dir

from typing import Optional

from tqdm import tqdm

from poliwag.antismash_processing.antismash_json import as_to_polymers
from poliwag.antismash_processing.polymer import PredictedPolymer


def parse_arguments() -> Namespace:
    parser = ArgumentParser(description='Get polymers from antismash output')
    parser.add_argument('-as', '--antismash_output', type=str,
                        help='Path to a folder containing antiSMASH .json output files directly '
                        '(one file per genome), not a folder of per-genome subfolders')
    parser.add_argument('--per_gene', action="store_true", default=False,
                        help="If given, output polymers per gene instead of per candidate cluster")
    parser.add_argument('-o', '--out_dir', type=str, help='Path to output folder')

    args = parser.parse_args()
    return args


def _format_domain_ids(domain_ids: Optional[list[Optional[tuple[str, str]]]], length: int, separator: str) -> str:
    """Format a polymer's per-monomer domain provenance as a single "Domains"
    column, positionally aligned with the "Polymer" column (same separator,
    one token per monomer). Each token is ``locus_tag,domain_index`` (the
    same identifiers used in an amplicon FASTA header, see
    ``amp_binding_domains_to_fasta.AmpBindingDomain.fasta_id``), or ``NA`` for
    a monomer with no corresponding AMP-binding domain (e.g. a PKS_AT/CAL
    building block), or when ``domain_ids`` wasn't recorded at all.
    """
    if not domain_ids:
        return separator.join(["NA"] * length)
    tokens = []
    for entry in domain_ids:
        if entry is None:
            tokens.append("NA")
        else:
            locus_tag, domain_index = entry
            tokens.append(f"{locus_tag},{domain_index}")
    return separator.join(tokens)


def _parse_domain_ids(domains_field: str, separator: str) -> list[Optional[tuple[str, str]]]:
    """Inverse of ``_format_domain_ids``."""
    domain_ids: list[Optional[tuple[str, str]]] = []
    for token in domains_field.strip().split(separator):
        if token == "NA" or not token:
            domain_ids.append(None)
        else:
            locus_tag, domain_index = token.split(",", 1)
            domain_ids.append((locus_tag, domain_index))
    return domain_ids


def write_polymers(polymers: list[PredictedPolymer], out_file: str, write_candidate_clusters: bool = True, separator: str = '|') -> None:
    """Write polymers to output file

    :param polymers: list of polymers
    :param out_file: path to output file
    :param write_candidate_clusters: whether to write the candidate cluster number for the polymer or not
    :param separator: separator to separate polymers

    Also writes a trailing "Domains" column with each polymer's per-monomer
    domain provenance (see ``_format_domain_ids``), so a polymer position can
    later be matched back to the exact AMP-binding domain amplicon that
    encodes it, instead of only to the candidate cluster/gene it came from.
    """

    with open(out_file, 'w') as out:
        if write_candidate_clusters:
            out.write("Accession\tRegion\tCandidate\tPolymer\tDomains\n")
        else:
            out.write("Accession\tRegion\tLocus tag\tPolymer\tDomains\n")
        for polymer in polymers:
            domains_column = _format_domain_ids(polymer.domain_ids, len(polymer.polymer), separator)
            if write_candidate_clusters:
                out.write(f"{polymer.accession}\t{polymer.region}\t{polymer.candidate_cluster}\t{separator.join(polymer.polymer)}\t{domains_column}\n")
            else:
                out.write(f"{polymer.accession}\t{polymer.region}\t{polymer.locus_tag}\t{separator.join(polymer.polymer)}\t{domains_column}\n")

def parse_polymers_bulk(antismash_json_folder: str, per_gene: bool = False) -> list[PredictedPolymer]:
    """Parse polymers from a flat folder of antiSMASH JSON output files (must
    be a version of antiSMASH with PARAS integrated).

    :param antismash_json_folder: path to a folder containing antiSMASH .json
        output files directly (one file per genome) -- not a folder of
        per-genome subfolders (matches poliwag.amplicons.amp_binding_domains_to_fasta's
        expected layout, so the same antiSMASH output folder can feed both)
    :param per_gene: whether to output polymers per gene instead of per candidate cluster
    """
    predicted_polymers = []

    json_files = list(iterate_over_dir(antismash_json_folder, ".json"))
    for genome_name, file_path in tqdm(json_files, desc="Parsing polymers", unit="file"):
        try:
            polymers = as_to_polymers(file_path, per_gene)
            predicted_polymers.extend(polymers)
        except Exception:
            print(f"Could not decode input file {genome_name}")
            print(traceback.format_exc())

    return predicted_polymers


def read_unique_polymers(polymers_file: str, separator: str = '|') -> list[list[str]]:
    """Read polymers from a file with one polymer per line

    :param polymers_file: path to polymers file
    :param separator: separator to separate polymers into monomers
    :return: list of polymers

    """
    unique_polymers: list[list[str]] = []
    with open(polymers_file, 'r') as polymers:
        for polymer in polymers:
            polymer = polymer.strip()
            monomers = polymer.split(separator)
            unique_polymers.append(monomers)
    return unique_polymers

def read_polymers(polymers_file: str, separator: str = '|') -> list[PredictedPolymer]:
    """Read polymers from a tabular file containing genome context information

    :param polymers_file: path to polymers file
    :param separator: separator to separate polymers into monomers
    :return: list of polymers

    Reads a trailing "Domains" column if the file has one (as written by
    ``write_polymers``); older files without it are read the same as before,
    with ``PredictedPolymer.domain_ids`` left as ``None``.
    """
    predicted_polymers = []
    with open(polymers_file, 'r') as polymers:
        header = polymers.readline()
        write_candidate_clusters = "Candidate" in header
        has_domains_column = "Domains" in header
        for line in polymers:
            fields = line.rstrip('\n').split('\t')
            if write_candidate_clusters:
                accession, region, candidate, polymer = fields[:4]
                locus_tag = None
            else:
                accession, region, locus_tag, polymer = fields[:4]
                candidate = None

            domain_ids = None
            if has_domains_column and len(fields) > 4:
                domain_ids = _parse_domain_ids(fields[4], separator)

            polymer = PredictedPolymer(
                polymer.strip().split(separator), accession, int(region), candidate, locus_tag, domain_ids
            )
            predicted_polymers.append(polymer)

    return predicted_polymers


def get_unique_polymers(predicted_polymers: list[PredictedPolymer]) -> list[str]:
    """
    Return a list of unique polymers from a list of predicted polymers

    :param predicted_polymers: list of predicted polymers
    :return: list of unique polymers
    """
    unique_polymers = set()
    for polymer in predicted_polymers:
        unique_polymers.add(str(polymer))

    unique_polymers = list(unique_polymers)
    unique_polymers.sort()

    return unique_polymers

def main() -> None:
    args = parse_arguments()
    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)

    polymers = parse_polymers_bulk(args.antismash_output,args.per_gene)
    out_all = os.path.join(args.out_dir, 'all_polymers.txt')
    out_unique = os.path.join(args.out_dir, 'unique_polymers.txt')
    if args.per_gene:
        write_candidate_clusters = False
    else:
        write_candidate_clusters = True

    write_polymers(polymers, out_all, write_candidate_clusters)

    with open(out_unique, 'w') as out:
        for polymer in get_unique_polymers(polymers):
            out.write(f"{polymer}\n")

if __name__ == '__main__':

    main()