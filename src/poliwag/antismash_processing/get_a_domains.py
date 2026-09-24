import sys
from argparse import ArgumentParser, Namespace
import traceback


from poliwag.antismash_processing.antismash_json import as_to_a_domains
from poliwag.antismash_processing.domain import AdenylationDomain
from poliwag.utils import iterate_over_dir


def parse_arguments() -> Namespace:
    parser = ArgumentParser(description='Get polymers from antismash output')
    parser.add_argument('-as', '--antismash_output', type=str, help='Path to antismash output folder')
    parser.add_argument('-o', '--out_file', type=str, help='Path to output file')

    args = parser.parse_args()
    return args

def get_a_domains(antismash_output_folder: str) -> list[AdenylationDomain]:
    """Parse adenylation domains from antismash output and write to file

    :param antismash_output_folder: path to antismash output folder
    """
    all_a_domains: list[AdenylationDomain] = []

    counter = 0

    for _, folder_path in iterate_over_dir(antismash_output_folder, get_dirs=True):
        counter += 1
        for genome_name, file_path in iterate_over_dir(folder_path, '.json'):
            try:
                a_domains = as_to_a_domains(file_path)
                all_a_domains.extend(a_domains)

            except Exception:
                print(f"Could not decode input file {genome_name}")
                print(traceback.format_exc())

        if counter % 100 == 0:
            print(f"Processed {counter} antismash files")

    return all_a_domains

def main():
    args = parse_arguments()
    a_domains = get_a_domains(args.antismash_output)
    with open(args.out_file, 'w') as f:
        for a_domain in a_domains:
            a_domain.write_to_fasta(f)

if __name__ == "__main__":
    main()

