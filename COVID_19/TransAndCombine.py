
import argparse
import os
import subprocess
import re
from tqdm import tqdm

PDB_TRANS = "~dina/utils/pdb_trans "
SPIKE_PDB = "/cs/usr/tomer.cohen13/lab/nanobodies/COVID_19/6vxx.pdb "

def get_trans_pdb(folder):
    pdb_name = os.path.basename(folder)
    os.chdir(folder)
    i = 0
    with open("best_100_trans", 'r') as trans_file:
        trans = trans_file.readline().strip()
        while trans:
            nanobody_trans = "nanobody_trans_" + str(i)+".pdb"
            subprocess.run(PDB_TRANS + trans + " < " + pdb_name + "_nanobody.pdb > " + nanobody_trans, shell=True)
            subprocess.run("cat " + SPIKE_PDB + " >> " + nanobody_trans, shell=True)
            trans = trans_file.readline().strip()
            i = i + 1
    os.chdir("..")


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="directory path containing the pdb directories")
    args = parser.parse_args()
    os.chdir(args.directory)

    for directory in tqdm(os.listdir(os.getcwd())):
        #  if the folder is pdb folder
        if os.path.isdir(directory) and re.fullmatch("[a-zA-Z0-9]{4}_[0-9]", directory):
            get_trans_pdb(directory)
