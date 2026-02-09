"""
Script to run antiSMASH on a large collection of genomes from the command line. Only region.gbk files are saved.
"""

import os
import subprocess
from typing import Optional
from shutil import move, rmtree
from argparse import ArgumentParser, Namespace

from poliwag.utils import iterate_over_dir


def parse_arguments() -> Namespace:
    """
    Parse command line arguments
    """
    parser = ArgumentParser(description="Run antiSMASH on downloaded ncbi dataset.")
    parser.add_argument('-i', '--input_directory', type=str, required=True,
                        help="Input directory. \
                        Folder architecture: one folder per genome contaiing a file called genomic.gbff")
    parser.add_argument('-o', '--output_directory', type=str, required=True,
                        help="Output directory. \
                        Will contain one folder per genome labelled by RefSeq assembly accession")
    parser.add_argument('-r', '--rerun', action="store_true",
                        help="If given, rerun analysis even if output folder with the same accession already exists.")
    parser.add_argument('-c', '--cpus', type=int, default=2,
                        help="Number of CPUs for running antiSMASH. A low number of CPUs is recommended on small\
                             machines to prevent memory contention.")
    args = parser.parse_args()

    return args


def run_antismash(input_file: str, output_dir: str, cpus: int = 2,
                  output_basename: Optional[str] = None) -> None:
    """
    Run antiSMASH from command line

    input_gbk: str, path to input .gbff file
    output_dir: str, path to output directory
    output_basename: str, prefix for file names within the output directory. Default: None

    """
    command = ['antismash', '--cpus', f'{cpus}', '--output-dir', output_dir, '--genefinding-tool',
               'prodigal', input_file]

    if output_basename is not None:
        command.extend(['--output-basename', output_basename])

    subprocess.call(command)


def move_antismash_region_files(input_folder: str, output_folder: str) -> None:
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    for file_name, file_path in iterate_over_dir(input_folder, '.gbk'):
        if '.region' in file_name:
            new_path = os.path.join(output_folder, file_name)
            move(file_path, new_path)


def main():
    """
    Run antiSMASH on NCBI dataset
    """
    args = parse_arguments()
    if not os.path.exists(args.output_directory):
        os.mkdir(args.output_directory)

    if not args.rerun:
        genomes_processed = sum(1 for entry in os.listdir(args.output_directory) if
                                os.path.isdir(os.path.join(args.output_directory, entry)))
    else:
        genomes_processed = 0

    temp_folder = os.path.join(args.output_directory, "temp")
    if not os.path.exists(temp_folder):
        os.mkdir(temp_folder)

    for accession, folder_path in iterate_over_dir(args.input_directory, get_dirs=True):
        genome_path = os.path.join(folder_path, "genomic.gbff")
        if os.path.exists(genome_path):
            output_folder = os.path.join(args.output_directory, accession)
            antismash_temp_output = os.path.join(temp_folder, accession)
            if not os.path.exists(output_folder) or args.rerun:
                run_antismash(genome_path, antismash_temp_output, args.cpus, accession)
                move_antismash_region_files(antismash_temp_output, output_folder)
                rmtree(antismash_temp_output)
                genomes_processed += 1
            else:
                print(f"AntiSMASH was already run on {accession}. Continuing..")
        else:
            print(f"No genomic.gbff file found for {accession}. Continuing..")

        if genomes_processed % 10 == 0 and genomes_processed != 0:
            print(f"Number of genomes processed: {genomes_processed}")

    rmtree(temp_folder)


if __name__ == "__main__":
    main()
