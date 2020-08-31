
import pandas as pd
import numpy as np
import os
import argparse
import subprocess


COL = ["type", "name", "dope_score", "soap_score", "rmsd"]
FULL_COL = ["type", "name", "dope_score", "soap_score", "rmsd", "length"]

LENGTH_PATH = "/cs/labs/dina/tomer.cohen13/lengths"
LENGTH_FILE = "nano_length.txt"

TOP_LOOP_N = 5
TOP_MODEL_N = 5

SCORE = "soap_score"


def get_scores_data(pdb_folder):
    """
    reads the scores.txt from the pdb_folder into df (column names = COL)
    :param pdb_folder: path of the pdb folder
    :return: df (len(COL) columns)
    """
    return pd.read_csv(os.path.join(pdb_folder, "scores.txt"), sep=" ", header=None, names=COL, usecols=[0,1,3,5,7])


def get_names_best_loops_models(pdb_folder):
    """

    """
    df = get_scores_data(pdb_folder)
    top_loop_names = pd.DataFrame.sort_values(df[df["type"] == "LOOP"], by=SCORE)[0:TOP_LOOP_N]["name"]
    top_model_names = pd.DataFrame.sort_values(df[df["type"] == "MODEL"], by=SCORE)[0:TOP_MODEL_N]["name"]

    return top_model_names.tolist(), top_loop_names.tolist()


def get_best_one_pdb(pdb_folder):

    models_names, loops_names = get_names_best_loops_models(pdb_folder)
    for i in range(TOP_LOOP_N):
        os.rename(os.path.join(pdb_folder, loops_names[i]), os.path.join(pdb_folder, "loop_" + str(i) + ".pdb"))
    for j in range(TOP_MODEL_N):
        os.rename(os.path.join(pdb_folder, models_names[j]), os.path.join(pdb_folder, "model_" + str(j)) + ".pdb")
    subprocess.run("tar -cf " + os.path.join(pdb_folder, "nanobodies.tar") + " " + os.path.join(pdb_folder, "NANO*"), shell=True)
    subprocess.run("rm -f " + os.path.join(pdb_folder, "NANO*"), shell=True)


def get_best_pdbs(directory):

    for pdb_folder in os.listdir(directory):
        if os.path.isdir(pdb_folder) and os.path.exists(os.path.join(pdb_folder, "scores.txt")):
            get_best_one_pdb(pdb_folder)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="dirctory path containing the pdb directories")
    args = parser.parse_args()
    os.chdir(args.directory)

    get_best_pdbs(args.directory)




