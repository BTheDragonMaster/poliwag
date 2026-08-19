from argparse import ArgumentParser, Namespace
from pathlib import Path

from antismash.modules.nrps_pks.name_mappings import get_substrate_by_name

from poliwag.antismash_processing.get_polymers import read_unique_polymers

def parse_args() -> Namespace:
    parser = ArgumentParser(description="Preprocess polymers for alignment with NATU")
    parser.add_argument('-i', '--input', required=True, type=str,
                        help="Path to input polymers file")
    # parser.add_argument('-s', '--substrates',
    #                     required=True,
    #                     type=Path,
    #                     help="Path to substrates SMILES file")

    parser.add_argument('-o', '--output', required=True, type=str,
                        help="Path to output polymers file")

    args = parser.parse_args()
    return args

def get_full_name(monomer: str) -> str:
    is_d = False
    has_nme = False
    has_me = False
    if monomer.startswith('D-') and 'D-Lya' not in monomer:
        name = monomer[2:]
        is_d = True
    else:
        name = monomer

    if name.startswith('NMe-'):
        has_nme = True
        name = name[4:]
    if name.startswith('Me-'):
        has_me = True
        name = name[3:]

    try:
        long_name = get_substrate_by_name(name).long
    except ValueError:
        print(monomer)
        long_name = name

    if has_me:
        long_name = long_name # NATU currently does not handle C-methylations
    if has_nme:
        long_name = f"NMe-{long_name}"
    if is_d:
        long_name = f"D-{long_name}"
    return long_name

def convert_to_natu_format(polymer: list[str]) -> list[str]:
    natu_formatted_polymer = []
    for monomer in polymer:
        natu_formatted_polymer.append(get_full_name(monomer))


    return natu_formatted_polymer

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