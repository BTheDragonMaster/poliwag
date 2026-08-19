from sys import argv


def polymers_to_fasta(polymer_in: str, fasta_out: str) -> None:
    """Convert file of unique polymers to FASTA format"""
    with open(fasta_out, 'w') as out:
        with open(polymer_in, 'r') as polymers:
            polymer_number = 1
            for line in polymers:
                line = line.strip()
                if line:
                    sequence = line.split(' | ')
                    seq_string = '|'.join(sequence)
                    out.write(f">polymer_{polymer_number}\n{seq_string}\n")
                    polymer_number += 1


if __name__ == "__main__":
    polymers_to_fasta(argv[1], argv[2])