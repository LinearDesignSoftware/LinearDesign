from concurrent.futures import ProcessPoolExecutor
import subprocess as sp
from multiprocessing import cpu_count
import pathlib
import re

DATAPATH = "./data/proteins/CRISPR_protein.fasta"
DESIGNPATH = "./designs/proteins"
# LAMBDA = [i for i in range(0, 100, 5)]
LAMBDA = [0]


def subprocess_lineardesign(cmd1, cmd2):
    # Create the pipeline by connecting the two commands using a pipe
    p1 = sp.Popen(cmd1, stdout=sp.PIPE)
    p2 = sp.Popen(cmd2, stdin=p1.stdout, stdout=sp.PIPE)

    # Run the pipeline and capture the output
    output, error = p2.communicate()

    return output


if __name__ == "__main__":
    path = pathlib.Path(DATAPATH)
    split_path = path.parent / "split"
    split_path.mkdir(parents=True, exist_ok=True)

    sp.run(["split", "-l", "2", DATAPATH, f"./{str(split_path)}/"])
    items = list(pathlib.Path(split_path).glob("*"))

    with ProcessPoolExecutor(max_workers=cpu_count()) as executor:
        # output redirection from stdout to file
        for lambda_ in LAMBDA:
            lambda_group = []
            for item in items:
                # Define the first command in the pipeline
                cmd1 = ["cat", f"{str(item)}"]
                # Define the second command in the pipeline
                cmd2 = ["./lineardesign", "-l", str(lambda_)]

                future = executor.submit(subprocess_lineardesign, cmd1, cmd2)
                lambda_group.append(future)

            pattern = "j=\d*\\r"
            lambda_group = [future.result().decode() for future in lambda_group]
            lambda_group = map(lambda x: re.sub(pattern, "", x), lambda_group)  # Remove iteration mark

            # # Print the output
            # print(output)

            pathlib.Path(DESIGNPATH).mkdir(parents=True, exist_ok=True)
            with open(f"{DESIGNPATH}/{path.name}+lambda_{lambda_}.txt", "w") as f:
                f.writelines(lambda_group)
