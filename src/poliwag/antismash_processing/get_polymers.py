"""
Get predicted NRPS and PKS polymers from antiSMASH output
"""

import os
from argparse import ArgumentParser, Namespace
import traceback

from poliwag.utils import iterate_over_dir

from poliwag.antismash_processing.json import as_to_polymers
from poliwag.antismash_processing.polymer import PredictedPolymer


def parse_arguments() -> Namespace:
    parser = ArgumentParser(description='Get polymers from antismash output')
    parser.add_argument('-as', '--antismash_output', type=str, help='Path to antismash output folder')
    parser.add_argument('--per_gene', action="store_true", default=False,
                        help="If given, output polymers per gene instead of per candidate cluster")
    parser.add_argument('-o', '--out_dir', type=str, help='Path to output folder')

    args = parser.parse_args()
    return args


def write_polymers(polymers: list[PredictedPolymer], out_file: str, write_candidate_clusters: bool = True) -> None:
    """Write polymers to output file

    :param polymers: list of polymers
    :param out_file: path to output file
    :param write_candidate_clusters: whether to write the candidate cluster number for the polymer or not

    """

    with open(out_file, 'w') as out:
        if write_candidate_clusters:
            out.write("Accession\tRegion\tCandidate\tPolymer\n")
        else:
            out.write("Accession\tRegion\tLocus tag\tPolymer\n")
        for polymer in polymers:
            if write_candidate_clusters:
                out.write(f"{polymer.accession}\t{polymer.region}\t{polymer.candidate_cluster}\t{polymer}\n")
            else:
                out.write(f"{polymer.accession}\t{polymer.region}\t{polymer.locus_tag}\t{polymer}\n")

def parse_polymers_bulk(antismash_output_folder: str, per_gene: bool = False) -> list[PredictedPolymer]:
    """Parse polymers from antismash output (must be a version of antiSMASH with PARAS integrated)

    :param antismash_output_folder: path to antismash output folder
    :param per_gene: whether to output polymers per gene instead of per candidate cluster
    """
    predicted_polymers = []

    counter = 0

    for _, folder_path in iterate_over_dir(antismash_output_folder, get_dirs=True):
        counter += 1
        for genome_name, file_path in iterate_over_dir(folder_path, '.json'):
            print(genome_name)
            try:
                polymers = as_to_polymers(file_path, per_gene)
                predicted_polymers.extend(polymers)

            except Exception:
                print(f"Could not decode input file {genome_name}")
                print(traceback.format_exc())

        if counter % 100 == 0:
            print(f"Processed {counter} antismash files")

    return predicted_polymers


def read_unique_polymers(polymers_file: str) -> list[list[str]]:
    """Read polymers from a file with one polymer per line

    :param polymers_file: path to polymers file
    :return: list of polymers

    """
    unique_polymers: list[list[str]] = []
    with open(polymers_file, 'r') as polymers:
        for polymer in polymers:
            polymer = polymer.strip()
            monomers = polymer.split(' | ')
            unique_polymers.append(monomers)
    return unique_polymers

def read_polymers(polymers_file: str) -> list[PredictedPolymer]:
    """Read polymers from a tabular file containing genome context information

    :param polymers_file: path to polymers file
    :return: list of polymers
    """
    predicted_polymers = []
    with open(polymers_file, 'r') as polymers:
        header = polymers.readline()
        if "Candidate" in header:
            write_candidate_clusters = True
        else:
            write_candidate_clusters = False
        for line in polymers:
            if write_candidate_clusters:
                accession, region, candidate, polymer = line.split('\t')
                locus_tag = None
            else:
                accession, region, locus_tag, polymer = line.split('\t')
                candidate = None
            polymer = PredictedPolymer(polymer.strip().split(' | '), accession, int(region), candidate, locus_tag)
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