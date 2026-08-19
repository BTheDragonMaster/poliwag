import os
from subprocess import check_call
from argparse import ArgumentParser, Namespace

from Bio import SearchIO
from Bio.SearchIO._model import HSP

STARTING_POSITION = 138


def parse_arguments() -> Namespace:
    parser = ArgumentParser(description='Run in silico amplicon sequencing on fasta file containing A domain sequences')
    parser.add_argument("-f", "--fasta", type=str, required=True,
                        help="Path to fasta file containing A domain sequences")
    parser.add_argument("-o", "--output", type=str, required=True, help="Path to output fasta file")
    parser.add_argument("-d", "--hmm_database", type=str, required=True, help="Path to hmm database")
    parser.add_argument("-t", "--temp_dir", type=str, default=os.path.join(os.getcwd(), "tmp"),
                        help="Path to temp directory")
    args = parser.parse_args()
    return args

def run_hmmscan(fasta_file: str, hmm_database_dir: str, out_file: str) -> None:
    command = ["hmmscan", "-o", out_file, hmm_database_dir, fasta_file]
    check_call(command)


def build_alignment_map(hsp):
    """
    Maps alignment columns to query residue indices.
    """

    hit_i = hsp.hit_start - 1
    query_i = -1
    aln_to_query = {}

    for aln_pos, q_res in enumerate(hsp.query.seq):
        h_res = hsp.hit.seq[aln_pos]
        if h_res != ".":
            hit_i += 1

        if q_res != "-":
            query_i += 1
            aln_to_query[hit_i] = query_i


    return aln_to_query

def parse_hmm_results(path_in: str, hmmer_version: int = 3) -> dict[str, HSP]:
    """Parse hmmpfam2 output file and return dictionary of domain identifier to Biopython HSP instance.

    :param path_in: path to hmmpfam2 output file (hmmer-2).
    :type path_in: str
    :param hmmer_version: version of HMMer
    "type hmmer_version: int, default: 2. Must be 2 or 3
    :return: Dictionary mapping domain identifier to Biopython HSP instance.
    :rtype: Dict[str, HSP]
    """

    if hmmer_version not in [2, 3]:
        raise ValueError(f"Unknown HMMer version: {hmmer_version}")
    filtered_hits = {}

    hmmer_string = f"hmmer{hmmer_version}-text"

    # parse relevant information from hmmpfam2 output
    for result in SearchIO.parse(path_in, hmmer_string):
        for hsp in result.hsps:

            # filter hits based on bitscore and hit_id
            if hsp.bitscore > 20:
                if hsp.hit_id == "AMP-binding" or hsp.hit_id == "AMP-binding_C":

                    header = f"{result.id}|{hsp.hit_id}|{hsp.query_start}-{hsp.query_end}"
                    filtered_hits[header] = hsp

    return filtered_hits


def extract_amplicon(hsp, query_pos, amplicon_size=33):

    seq = hsp.query.seq.replace("-", "")

    return seq[query_pos:query_pos + amplicon_size]


def main():
    args = parse_arguments()
    if not os.path.exists(args.temp_dir):
        os.mkdir(args.temp_dir)

    tmp_hmm = os.path.join(args.temp_dir, "tmp.hmm")
    run_hmmscan(args.fasta, args.hmm_database, tmp_hmm)
    id_to_hit = parse_hmm_results(tmp_hmm)
    counter = 0
    with open(args.output, 'w') as out:
        for seq_id, hit in id_to_hit.items():
            aln_to_query = build_alignment_map(hit)
            if STARTING_POSITION in aln_to_query:
                amplicon = extract_amplicon(hit, aln_to_query[STARTING_POSITION], amplicon_size=33)
            else:
                print(f"Could not extract amplicon for sequence {seq_id}")
                counter += 1

            print(amplicon)
            out.write(f">{seq_id}\n{amplicon}\n")

    print(f"Could not find amplicons for {counter} sequences.")

if __name__ == "__main__":
    main()





