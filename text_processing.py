# Description: This file contains the functions for text processing


def input_preprocessing(path: str, mode=None) -> str:
    """
    Preprocess the input file to remove the unnecessary line breaks
    """
    with open(path, "r") as f:
        lines = f.read().split(">")[1:]  # Remove the first empty string
        buffer = []
        for line in lines:
            temp = line.split("\n")
            if mode is None:
                buffer.append(f'>{temp[0]}\n{"".join(temp[1:])}*\n')
            elif mode == "RNA":
                buffer.append(f'>{temp[0]}\n{"".join(temp[1:])}\n')

        preprocessed_file_path = f"{path}_preprocessed.fasta"
        with open(preprocessed_file_path, "w") as d:
            d.writelines(buffer)
        return preprocessed_file_path
