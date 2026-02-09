import os


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
