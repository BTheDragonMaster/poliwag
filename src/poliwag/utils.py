import os


def parse_fasta_file(path_in: str, separator: str = "") -> dict[str, str]:
    """Read fasta file and parses it into a dictionary format.

    :param path_in: Path to fasta file.
    :param separator: Separator to separate sequence elements
    :return: Dictionary of fasta sequences, where the key is the sequence header
        and the value is the sequence.
    :rtype: Dict[str, str]
    :raises FileNotFoundError: If the file at the specified path does not exist.
    """
    # check if the file exists
    if not os.path.exists(path_in):
        raise FileNotFoundError(f"File not found: {path_in}")

    # initialize the dictionary
    fasta_dict: dict[str, str] = {}

    with open(path_in, "r") as fo:

        # initialize the header and sequence list
        header = ""
        sequence: list[str] = []

        for line in fo:
            line = line.strip()

            if line.startswith(">"):
                if sequence:
                    # if the sequence list is not empty, join the list into a string
                    # and add it to the dictionary with the header as the key
                    fasta_dict[header] = separator.join(sequence)

                    header = line[1:]  # remove the ">" from the header
                    sequence = []  # reset the sequence list
                else:
                    header = line[1:]

            else:
                sequence.append(line)

        # if the header is not None, this means there is a sequence that has not
        # yet been added to the dictionary
        if len(header) > 0:
            fasta_dict[header] = separator.join(sequence)

    return fasta_dict

def iterate_over_dir(directory, suffix=None, get_dirs=False):
    for file_name in os.listdir(directory):
        path = os.path.join(directory, file_name)
        if suffix and file_name.endswith(suffix):
            if os.path.isfile(path) and not get_dirs:
                yield file_name, path
            elif os.path.isdir(path) and get_dirs:
                yield file_name, path
        elif suffix is None:
            if os.path.isfile(path) and not get_dirs:
                yield file_name, path
            elif os.path.isdir(path) and get_dirs:
                yield file_name, path
