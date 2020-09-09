
import os
import argparse
import re
import subprocess
from tqdm import tqdm

CLUSTER = "/cs/labs/dina/tomer.cohen13/InterfaceClustering/interface_cluster.linux "


def cluster(directory):
    os.chdir(directory)
    subprocess.run(CLUSTER + " -f antigen.pdb ref.pdb dock_data.csv 4 soap_score_cluster.res", shell=True)
    os.chdir("..")


if __name__ == '__main__':

    """
    runs the dock_pdb() function on all the pdbs folders in the given directory.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="directory path containing the pdb directories")
    args = parser.parse_args()
    os.chdir(args.directory)

    for directory in tqdm(os.listdir(args.directory)):
        #  if the folder is pdb folder
        if os.path.isdir(directory) and re.fullmatch("[a-zA-Z0-9]{4}_[0-9]", directory):
            cluster(os.path.join(args.directory, directory))

