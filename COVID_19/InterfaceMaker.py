import argparse
import os
import subprocess
import re
import pandas as pd

INTERFACE = "/cs/labs/dina/tomer.cohen13/nanobodies/epiDock/interface "
HEADER = "/cs/usr/tomer.cohen13/lab/nanobodies/COVID_19/header.csv"


def get_interface(folder):
    os.chdir(folder)
    subprocess.run(INTERFACE + " ABC H 6 nanobody_trans_*.pdb", shell=True)
    df = pd.read_csv("epi.csv", header=None)
    os.chdir("..")
    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="directory path containing the pdb directories")
    args = parser.parse_args()
    os.chdir(args.directory)

    with open("epiDock_summey.csv", "w") as epi_file:
        header_df = pd.read_csv(HEADER)
        header_df.insert(0, column="PDB", value=None)
        header_df.to_csv(epi_file, header=True, index=False)
        for directory in os.listdir(os.getcwd()):
            #  if the folder is pdb folder
            if os.path.isdir(directory) and re.fullmatch("[a-zA-Z0-9]{4}_[0-9]", directory):
                epi_df = get_interface(directory)
                epi_df.insert(0, column="PDB", value=directory)
                epi_df.to_csv(epi_file, header=False, index=False)