from argparse import ArgumentParser, Namespace

from poliwag.alignment.natu import convert_to_natu_format

from poliwag.antismash_processing.get_polymers import read_polymers, write_polymers

def parse_args() -> Namespace:
    parser = ArgumentParser(description="Preprocess polymers for alignment with NATU")
    parser.add_argument('-i', '--input', required=True, type=str,
                        help="Path to input polymers file")
    parser.add_argument('-o', '--output', required=True, type=str,
                        help="Path to output polymers file")
    parser.add_argument('-g', '--per_gene', action='store_true',
                        help="If given, input polymer file contains per-gene polymers")

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    polymers = read_polymers(args.input)
    for polymer in polymers:
        natu_polymer = convert_to_natu_format(polymer.polymer)
        polymer.polymer = natu_polymer

    if args.per_gene:
        write_polymers(polymers, args.output, False, separator='|')
    else:
        write_polymers(polymers, args.output, True, separator='|')

if __name__ == "__main__":
    main()