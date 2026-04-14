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
    parser.add_argument('-i', '--input', type=str, required=True,
                        help="Input directory. \
                        Folder architecture: one folder per genome containing a file called genomic.gbff")
    parser.add_argument('--from_ncbi_dataset', action='store_true',
                        help="If given, the input is an NCBI fetch.txt file, obtained by running the fetch_genomes script.\
                             Genomes will be downloaded 100 at a time and erased after running antiSMASH.")
    parser.add_argument('-o', '--output_directory', type=str, required=True,
                        help="Output directory. \
                        Will contain one folder per genome labelled by RefSeq assembly accession")
    parser.add_argument('-r', '--rerun', action="store_true",
                        help="If given, rerun analysis even if output folder with the same accession already exists.")
    parser.add_argument('-c', '--cpus', type=int, default=2,
                        help="Number of CPUs for running antiSMASH. A low number of CPUs is recommended on small\
                             machines to prevent memory contention.")
    parser.add_argument('-js', )
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

def move_antismash_json_file(input_folder: str, output_folder: str) -> None:
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    for file_name, file_path in iterate_over_dir(input_folder, '.json'):
        new_path = os.path.join(output_folder, file_name)
        move(file_path, new_path)

def get_genomes_processed(output_dir: str, rerun: bool) -> int:
    if not rerun:
        genomes_processed = sum(1 for entry in os.listdir(output_dir) if
                                os.path.isdir(os.path.join(output_dir, entry)))
    else:
        genomes_processed = 0

    return genomes_processed

def run_antismash_from_ncbi_query(fetch_ncbi_file: str, temp_folder: str, output_folder: str,
                                  rerun: bool, cpus: int) -> None:
    genomes_processed = get_genomes_processed(output_folder, rerun)
    raise NotImplementedError


def run_antismash_from_genomes_folder(input_folder: str, temp_folder: str, antismash_output_folder: str, rerun: bool,
                                      cpus: int) -> None:
    genomes_processed = get_genomes_processed(antismash_output_folder, rerun)
    for accession, folder_path in iterate_over_dir(input_folder, get_dirs=True):
        genome_path = os.path.join(folder_path, "genomic.gbff")
        if os.path.exists(genome_path):
            output_folder = os.path.join(antismash_output_folder, accession)
            antismash_temp_output = os.path.join(temp_folder, accession)
            if not os.path.exists(output_folder) or rerun:
                run_antismash(genome_path, antismash_temp_output, cpus, accession)
                # move_antismash_region_files(antismash_temp_output, output_folder)
                move_antismash_json_file(antismash_temp_output, output_folder)
                rmtree(antismash_temp_output)
                genomes_processed += 1
            else:
                print(f"AntiSMASH was already run on {accession}. Continuing..")
        else:
            print(f"No genomic.gbff file found for {accession}. Continuing..")

        if genomes_processed % 10 == 0 and genomes_processed != 0:
            print(f"Number of genomes processed: {genomes_processed}")


def main():
    """
    Run antiSMASH on NCBI dataset
    """
    args = parse_arguments()
    if not os.path.exists(args.output_directory):
        os.mkdir(args.output_directory)

    temp_folder = os.path.join(args.output_directory, "temp")
    if not os.path.exists(temp_folder):
        os.mkdir(temp_folder)

    if not args.from_ncbi_dataset:
        run_antismash_from_genomes_folder(args.input, temp_folder, args.output_directory, args.rerun, args.cpus)
    else:
        run_antismash_from_ncbi_query(args.input, temp_folder, args.output_directory, args.rerun, args.cpus)

    rmtree(temp_folder)


if __name__ == "__main__":
    main()
