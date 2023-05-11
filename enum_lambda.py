from concurrent.futures import ProcessPoolExecutor
import subprocess as sp
from multiprocessing import cpu_count
import pathlib
import re
import os

DATAPATH = "./data/proteins/testseq"
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


def split_directory_cleanup(path: pathlib.Path) -> None:
    # get a list of all files in the directory
    file_list = os.listdir(str(path))

    # iterate over each file in the directory and remove it
    for filename in file_list:
        file_path = os.path.join(str(path), filename)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
                # print(f"Deleted {file_path}")
        except Exception as e:
            print(f"Failed to delete {file_path}. Reason: {e}")


if __name__ == "__main__":
    path = pathlib.Path(DATAPATH)
    split_path = path.parent / "split"
    try:
        split_path.mkdir(parents=True)
    except FileExistsError:
        print(f"Directory {split_path} already exists.")
        print(f"Cleaning up {split_path}...")
        split_directory_cleanup(split_path)

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
            lambda_group = map(
                lambda x: re.sub(pattern, "", x), lambda_group
            )  # Remove iteration mark

            # # Print the output
            # print(output)

            pathlib.Path(DESIGNPATH).mkdir(parents=True, exist_ok=True)
            with open(f"{DESIGNPATH}/{path.name}+lambda_{lambda_}.txt", "w") as f:
                f.writelines(lambda_group)
