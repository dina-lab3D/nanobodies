
import os
import argparse
import re
import subprocess
from tqdm import tqdm


# interface_cluster path
CLUSTER = "/cs/labs/dina/tomer.cohen13/InterfaceClustering/interface_cluster.linux "


def cluster(folder):
    """
    runs interface_cluster on the antigen.pdb ref.pdb dock_data.csv file in folder, put the results in soap_score_cluster.res
    :param folder: the pdb folder (with dock_data.csv file)
    :return: None
    """
    os.chdir(folder)
    subprocess.run(CLUSTER + " -f antigen.pdb ref.pdb dock_data.csv 4 soap_score_cluster.res", shell=True)
    os.chdir("..")


if __name__ == '__main__':

    """
    runs interface_cluster.linux on all the pdb folders in the directory (must have a dock_data.csv file from PyDock)
    saves the score results after clustering in soap_score_cluster.res
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="directory path containing the pdb directories")
    args = parser.parse_args()
    os.chdir(args.directory)

    for directory in tqdm(os.listdir(args.directory)):
        #  if the folder is pdb folder
        if os.path.isdir(directory) and re.fullmatch("[a-zA-Z0-9]{4}_[0-9]", directory):
            cluster(os.path.join(args.directory, directory))

