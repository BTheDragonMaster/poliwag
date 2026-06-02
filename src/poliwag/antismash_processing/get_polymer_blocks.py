import os
from argparse import ArgumentParser, Namespace
import traceback

from poliwag.utils import iterate_over_dir

from poliwag.antismash_processing.parse_json import polymers_from_json, parse_json
from poliwag.antismash_processing.polymer import PredictedPolymer


def parse_arguments() -> Namespace:
    parser = ArgumentParser(description='Get polymers from antismash output')
    parser.add_argument('-as', '--antismash_output', type=str, help='Path to antismash output folder')
    parser.add_argument('-o', '--out_dir', type=str, help='Path to output folder')

    args = parser.parse_args()
    return args


def write_polymers_from_folder(antismash_output_folder: str, out_file: str) -> None:
    with open(out_file, 'w') as out:
        for _, folder_path in iterate_over_dir(antismash_output_folder, get_dirs=True):
            for genome_name, file_path in iterate_over_dir(folder_path, '.json'):
                try:
                    polymers = polymers_from_json(file_path)
                    for polymer in polymers:
                        out.write(f"{polymer.accession}\t{polymer.region}\t{polymer.candidate_cluster}\t{polymer}\n")
                except Exception as e:
                    print(f"Could not decode input file {genome_name}")
                    print(e)

def write_polymers(polymers: list[PredictedPolymer], out_file: str) -> None:
    with open(out_file, 'w') as out:
        for polymer in polymers:
            out.write(f"{polymer.accession}\t{polymer.region}\t{polymer.candidate_cluster}\t{polymer}\n")

def parse_polymers_bulk(antismash_output_folder: str) -> list[PredictedPolymer]:
    predicted_polymers = []

    counter = 0

    for _, folder_path in iterate_over_dir(antismash_output_folder, get_dirs=True):
        counter += 1
        for genome_name, file_path in iterate_over_dir(folder_path, '.json'):
            print(genome_name)
            try:
                polymers = parse_json(file_path)
                predicted_polymers.extend(polymers)
            except Exception as e:
                print(f"Could not decode input file {genome_name}")
                print(traceback.format_exc())

        if counter % 100 == 0:
            print(f"Processed {counter} antismash files")

    return predicted_polymers


def read_polymers(polymers_file: str) -> list[PredictedPolymer]:
    predicted_polymers = []
    with open(polymers_file, 'r') as polymers:
        for line in polymers:
            accession, region, candidate, polymer = line.split('\t')
            polymer = PredictedPolymer(polymer.strip().split(' | '), accession, int(region), int(candidate))
            predicted_polymers.append(polymer)

    return predicted_polymers


def get_unique_polymers(predicted_polymers: list[PredictedPolymer]) -> list[str]:
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

    polymers = parse_polymers_bulk(args.antismash_output)
    out_all = os.path.join(args.out_dir, 'all_polymers.txt')
    out_unique = os.path.join(args.out_dir, 'unique_polymers.txt')
    write_polymers(polymers, out_all)

    with open(out_unique, 'w') as out:
        for polymer in get_unique_polymers(polymers):
            out.write(f"{polymer}\n")

if __name__ == '__main__':

    main()