from argparse import ArgumentParser, Namespace
import subprocess


def parse_arguments() -> Namespace:
    """
    Parse command line arguments
    """
    parser = ArgumentParser(description="Run antiSMASH on downloaded ncbi dataset.")
    parser.add_argument('-i', '--input_file', type=int, default=201174,
                        help="Taxon identifier of the genomes to analyse")
    parser.add_argument('-o', '--out_file', type=str, required=True,
                        help="Output file")

    args = parser.parse_args()

    return args


def fetch_ncbi_entries(taxon: int = 201174, out_file: str = "actinobacterial_refseq_gbffs.zip") -> None:
    """
    Run the ncbi-datasets-cli from the command line to fetch all genomes from the taxon of interest.

    taxon: int, default: 201174, actinomycetota
    out_file: str, default: actinobacterial_refseq_gbffs.zip
    """
    command = ['datasets', 'download', 'genome', 'taxon', f'{taxon}', '--assembly-source', 'refseq',
               '--include', 'gbff', '--dehydrated', '--filename', out_file]

    subprocess.call(command)

    command = ['unzip', out_file]

    subprocess.call(command)


def main():
    args = parse_arguments()

    fetch_ncbi_entries(args.taxon, args.out_file)


if __name__ == "__main__":
    main()
