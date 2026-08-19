from argparse import ArgumentParser, Namespace
from pathlib import Path

from natu.cli import read_monomer_fasta
import pandas as pd

def parse_args() -> Namespace:
    parser = ArgumentParser(description="Preprocess polymers for alignment with NATU")
    parser.add_argument('-i', '--input',
                        required=True,
                        type=Path,
                        help="Path to input polymers fasta")
    parser.add_argument('-n', '--natu_search_results',
                        required=True,
                        type=Path,
                        help="Path to natu search results")
    parser.add_argument('-o', '--output', required=True, type=Path,
                        help="Path to output polymers fasta")

    args = parser.parse_args()
    return args


def parse_natu_search_results(natu_search_results: Path) -> pd.DataFrame:

    df = pd.read_csv(natu_search_results, sep='\t')

    return df

def main():
    args = parse_args()
    search_results = parse_natu_search_results(args.natu_search_results)
    subjects = list(search_results["subject"])
    headers, sequences = zip(*read_monomer_fasta(args.input))

    with open(args.output, 'w') as out:

        for i, header in enumerate(headers):
            if header in subjects:
                out.write(f">{header}\n{'|'.join(sequences[i])}\n")

if __name__ == "__main__":
    main()
