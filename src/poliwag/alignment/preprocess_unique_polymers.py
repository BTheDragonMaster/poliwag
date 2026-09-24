from argparse import ArgumentParser, Namespace
from pathlib import Path

from poliwag.alignment.natu import convert_to_natu_format

from poliwag.antismash_processing.get_polymers import read_unique_polymers

def parse_args() -> Namespace:
    parser = ArgumentParser(description="Preprocess polymers for alignment with NATU")
    parser.add_argument('-i', '--input', required=True, type=str,
                        help="Path to input polymers file")

    parser.add_argument('-o', '--output', required=True, type=str,
                        help="Path to output polymers file")

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    unique_polymers = read_unique_polymers(args.input)
    natu_polymers: list[list[str]] = []
    natu_monomers: set[str] = set()
    for polymer in unique_polymers:
        natu_polymer = convert_to_natu_format(polymer)
        for monomer in natu_polymer:
            natu_monomers.add(monomer)
        natu_polymers.append(natu_polymer)

    with open(args.output, 'w') as out:
        for polymer in natu_polymers:
            out.write('|'.join(polymer))
            out.write('\n')

    for monomer in natu_monomers:
        if monomer.startswith('Me-'):
            print(monomer)


if __name__ == "__main__":
    main()