
import pandas as pd
import os
import argparse
import subprocess
from tqdm import tqdm


# columns names for the data frame
COL = ["type", "name", "dope_score", "soap_score", "rmsd", "cdr1_rmsd", "cdr2_rmsd", "cdr3_rmsd"]
FULL_COL = ["type", "name", "dope_score", "soap_score", "rmsd", "cdr1_rmsd", "cdr2_rmsd", "cdr3_rmsd", "length"]

# path to the directory containing all the lengths of the pdbs
LENGTH_PATH = "/cs/labs/dina/tomer.cohen13/lengths"

# length file name
LENGTH_FILE = "nano_length.txt"

# number of models to take from loops
TOP_LOOP_N = 5

# number of models to take from models
TOP_MODEL_N = 5

# score to use for choosing the best models
SCORE = "soap_score"

# for NanoNet training change to True
NANO_NET = True
NANO_NUM = 10



def get_scores_data(pdb_folder):
    """
    reads the scores.txt from the pdb_folder into df (column names = COL)
    :param pdb_folder: path of the pdb folder
    :return: df (len(COL) columns)
    """
    return pd.read_csv(os.path.join(pdb_folder, "scores.txt"), sep=" ", header=None, names=COL, usecols=[0,1,3,5,7,9,11,13])


def get_names_best_loops_models(pdb_folder):
    """
    returns 2 lists- first with top models names, second with top loops names
    :param pdb_folder: the pdb folder to use
    :return: 2 lists ( with length = TOP_MODEL_N, TOP_LOOP_N)
    """
    df = get_scores_data(pdb_folder)
    if NANO_NET:
        top_names = pd.DataFrame.sort_values(df, by="cdr1_rmsd")[0:NANO_NUM]["name"]
        # top_names["name"] = ["model_" + str(i) for i in range(NANO_NUM)]
        top_names.to_csv(os.path.join(pdb_folder, "top_models_rmsd.csv"))
        return top_names.tolist(),[]
    top_loop_names = pd.DataFrame.sort_values(df[df["type"] == "LOOP"], by=SCORE)[0:TOP_LOOP_N]["name"]
    top_model_names = pd.DataFrame.sort_values(df[df["type"] == "MODEL"], by=SCORE)[0:TOP_MODEL_N]["name"]

    return top_model_names.tolist(), top_loop_names.tolist()


def get_best_one_pdb(pdb_folder):
    """
    finds the best models generated for each type (model, loop) according to the score.
    renames them and sends the rest to tar file.
    :param pdb_folder: the pdb to use
    :return: None
    """
    models_names, loops_names = get_names_best_loops_models(pdb_folder)
    for i in range(len(loops_names)):
        os.rename(os.path.join(pdb_folder, loops_names[i]), os.path.join(pdb_folder, "loop_" + str(i) + ".pdb"))
    for j in range(len(models_names)):
        os.rename(os.path.join(pdb_folder, models_names[j]), os.path.join(pdb_folder, "model_" + str(j)) + ".pdb")
    subprocess.run("tar -cf " + os.path.join(pdb_folder, "nanobodies.tar") + " " + os.path.join(pdb_folder, "NANO*"), shell=True)
    subprocess.run("rm -f " + os.path.join(pdb_folder, "NANO*"), shell=True)


def get_best_pdbs(directory):
    """
    runs get_best_one_pdb on each pdb folder in directory
    :param directory: directory containing pdb folders
    :return: None
    """
    for pdb_folder in tqdm(os.listdir(directory)):
        if os.path.isdir(pdb_folder):
            if os.path.exists(os.path.join(pdb_folder, "scores.txt")) and not os.path.exists(os.path.join(pdb_folder, "nanobodies.tar")):
                get_best_one_pdb(pdb_folder)


if __name__ == '__main__':
    """
    this program chooses the top models generated by NanoModelScript.py according to the chosen score, changes their 
    names and sends the rest of the models to a tar file
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="dirctory path containing the pdb directories")
    args = parser.parse_args()
    os.chdir(args.directory)

    get_best_pdbs(args.directory)




