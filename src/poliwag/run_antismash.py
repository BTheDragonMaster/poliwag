"""
Script to run antiSMASH on a large collection of genomes from the command line. Only region.gbk files are saved.
"""

import os
import shlex
import subprocess
from typing import List, Optional
from shutil import move, rmtree, which
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
    parser.add_argument('--from_fna_folder', action='store_true',
                        help="If given, the input is a flat folder of .fna (nucleotide FASTA) files, one per genome. \
                             Each file is run through antiSMASH and results are written to a folder named after the file.")
    parser.add_argument('-o', '--output_directory', type=str, required=True,
                        help="Output directory. \
                        Will contain one folder per genome labelled by RefSeq assembly accession")
    parser.add_argument('-r', '--rerun', action="store_true",
                        help="If given, rerun analysis even if output folder with the same accession already exists.")
    parser.add_argument('-c', '--cpus', type=int, default=2,
                        help="Number of CPUs for running antiSMASH. A low number of CPUs is recommended on small\
                             machines to prevent memory contention.")
    parser.add_argument('-j', '--jobs', type=int, default=1,
                        help="Number of antiSMASH runs to execute at the same time using GNU parallel. Only available \
                             with --from_fna_folder. Total CPU use is roughly jobs x cpus. Default: 1 (sequential, \
                             GNU parallel not needed).")
    parser.add_argument('-js', )
    args = parser.parse_args()

    return args


def build_antismash_command(input_file: str, output_dir: str, cpus: int = 2,
                            output_basename: Optional[str] = None) -> List[str]:
    """
    Build the antiSMASH command line as a list of arguments
    """
    command = ['antismash', '--cpus', f'{cpus}', '--output-dir', output_dir, '--genefinding-tool',
               'prodigal', input_file]

    if output_basename is not None:
        command.extend(['--output-basename', output_basename])

    return command


def run_antismash(input_file: str, output_dir: str, cpus: int = 2,
                  output_basename: Optional[str] = None) -> int:
    """
    Run antiSMASH from command line. Returns the antiSMASH exit code (0 on success).

    input_gbk: str, path to input .gbff file
    output_dir: str, path to output directory
    output_basename: str, prefix for file names within the output directory. Default: None

    """
    return subprocess.call(build_antismash_command(input_file, output_dir, cpus, output_basename))


def run_commands_with_gnu_parallel(commands: List[List[str]], jobs: int, cpus: int, joblog_path: str) -> List[int]:
    """
    Run commands with GNU parallel, at most `jobs` at a time. Returns one exit code per command, in input order.
    A command killed by a signal gets the negative signal number, like subprocess does.

    Thread-count variables for math libraries are capped at `cpus` so that each job stays within its CPU budget.
    """
    if not commands:
        return []

    env = dict(os.environ)
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        env[variable] = str(cpus)

    if os.path.exists(joblog_path):
        os.remove(joblog_path)

    command_lines = "\n".join(shlex.join(command) for command in commands) + "\n"
    subprocess.run(['parallel', '--jobs', str(jobs), '--joblog', joblog_path],
                   input=command_lines, text=True, env=env)

    # joblog columns: Seq Host Starttime JobRuntime Send Receive Exitval Signal Command. Seq starts at 1
    return_codes = [-1] * len(commands)
    if os.path.exists(joblog_path):
        with open(joblog_path) as joblog:
            next(joblog, None)
            for line in joblog:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8 or not fields[0].isdigit():
                    continue
                index = int(fields[0]) - 1
                exit_value, signal = int(fields[6]), int(fields[7])
                if 0 <= index < len(commands):
                    return_codes[index] = -signal if signal else exit_value

    return return_codes


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


def collect_antismash_result(accession: str, return_code: int, antismash_temp_output: str,
                             output_folder: str) -> bool:
    """
    Move the antiSMASH json to its final folder and clean up the temp folder. If antiSMASH failed,
    report it and clean up instead. Returns True on success.
    """
    if return_code != 0 or not os.path.isdir(antismash_temp_output):
        print(f"AntiSMASH failed on {accession} (exit code {return_code}). Continuing..")
        if os.path.isdir(antismash_temp_output):
            rmtree(antismash_temp_output)
        return False

    move_antismash_json_file(antismash_temp_output, output_folder)
    rmtree(antismash_temp_output)
    return True


def run_antismash_from_fna_folder(input_folder: str, temp_folder: str, antismash_output_folder: str, rerun: bool,
                                  cpus: int, jobs: int = 1) -> None:
    """
    Run antiSMASH on every .fna file in a folder. Output for genome.fna goes in <antismash_output_folder>/genome

    With jobs > 1, up to `jobs` antiSMASH runs are executed at the same time with GNU parallel, each using `cpus` CPUs.
    """
    genomes_processed = get_genomes_processed(antismash_output_folder, rerun)

    pending = []
    for file_name, file_path in sorted(iterate_over_dir(input_folder, '.fna')):
        accession = file_name[:-len('.fna')]
        if not os.path.exists(os.path.join(antismash_output_folder, accession)) or rerun:
            pending.append((accession, file_path))
        else:
            print(f"AntiSMASH was already run on {accession}. Continuing..")

    if jobs > 1:
        print(f"Running {len(pending)} genomes, up to {jobs} at a time with {cpus} CPUs each "
              f"({jobs * cpus} CPUs in total).")
        commands = [build_antismash_command(file_path, os.path.join(temp_folder, accession), cpus, accession)
                    for accession, file_path in pending]
        joblog_path = os.path.join(temp_folder, "parallel_joblog.tsv")
        return_codes = run_commands_with_gnu_parallel(commands, jobs, cpus, joblog_path)
        succeeded = 0
        for (accession, _), return_code in zip(pending, return_codes):
            if collect_antismash_result(accession, return_code, os.path.join(temp_folder, accession),
                                        os.path.join(antismash_output_folder, accession)):
                succeeded += 1
        print(f"Finished: {succeeded} genomes succeeded, {len(pending) - succeeded} failed.")
        return

    for accession, file_path in pending:
        antismash_temp_output = os.path.join(temp_folder, accession)
        return_code = run_antismash(file_path, antismash_temp_output, cpus, accession)
        if collect_antismash_result(accession, return_code, antismash_temp_output,
                                    os.path.join(antismash_output_folder, accession)):
            genomes_processed += 1

        if genomes_processed % 10 == 0 and genomes_processed != 0:
            print(f"Number of genomes processed: {genomes_processed}")


def main():
    """
    Run antiSMASH on NCBI dataset
    """
    args = parse_arguments()

    if args.jobs < 1 or args.cpus < 1:
        raise SystemExit("--jobs and --cpus must be at least 1.")
    if args.jobs > 1 and not args.from_fna_folder:
        raise SystemExit("--jobs greater than 1 is only available with --from_fna_folder.")
    if args.jobs > 1 and which('parallel') is None:
        raise SystemExit("--jobs greater than 1 requires GNU parallel, but 'parallel' was not found on PATH.")

    if not os.path.exists(args.output_directory):
        os.mkdir(args.output_directory)

    temp_folder = os.path.join(args.output_directory, "temp")
    if not os.path.exists(temp_folder):
        os.mkdir(temp_folder)

    if args.from_fna_folder and args.from_ncbi_dataset:
        raise SystemExit("--from_fna_folder and --from_ncbi_dataset cannot be used together.")

    if args.from_fna_folder:
        run_antismash_from_fna_folder(args.input, temp_folder, args.output_directory, args.rerun, args.cpus,
                                      args.jobs)
    elif not args.from_ncbi_dataset:
        run_antismash_from_genomes_folder(args.input, temp_folder, args.output_directory, args.rerun, args.cpus)
    else:
        run_antismash_from_ncbi_query(args.input, temp_folder, args.output_directory, args.rerun, args.cpus)

    rmtree(temp_folder)


if __name__ == "__main__":
    main()
