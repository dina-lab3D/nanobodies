import pandas as pd
import numpy as np
import os, sys
import enum
import matplotlib.pyplot as plt


COL = ["type", "name", "2", "dope_score", "4", "soap_score","6" , "rmsd"]
FULL_COL = ["type", "name", "2", "dope_score", "4", "soap_score","6" , "rmsd", "length"]
LENGTH_PATH = "/cs/labs/dina/tomer.cohen13/lengths/"
LENGTH_FILE = "/nano_length.txt"


# Using enum class create enumerations
class ColInfo(enum.Enum):
    NAME = 1
    D_SCORE = 3
    S_SCORE = 5
    RMSD = 7
    LENGTH = 8


def generate_rmsd_graph(data, name):

    data["rmsd"] = pd.to_numeric(data["rmsd"])
    data["length"] = pd.to_numeric(data["length"])
    axes_boxplot = data.boxplot(column=["rmsd"], by=["length"], grid=False, showfliers=False)
    axes_boxplot.set_title("rmsd VS. length for " + name + " results")
    plt.savefig("rmsd_boxplot_" + name)



def generate_final_scores(folder_path):

    models, loops= [], []
    # create a list of the scores from the pdb files, divided into base model and sampled loops
    for file in os.listdir(folder_path):
        if os.path.isdir(folder_path + "/" + file) and os.path.exists(folder_path + "/" + file + "/scores.txt"):
            model, loop = parse_one_pdb(folder_path, file)
            models.append(model.tolist()[0])
            loops.append(loop.tolist()[0])
    return pd.DataFrame(data=models, columns=FULL_COL), pd.DataFrame(data=loops, columns=FULL_COL)


def parse_one_pdb(pdb_path, pdb_name):
    df = pd.read_csv(pdb_path + "/" + pdb_name + "/scores.txt", sep=" ", header=None, names=COL)

    loop_min_idx = df[df['type'] == "LOOP"]["rmsd"].idxmin()
    model_min_idx = df[df['type'] == "MODEL"]["rmsd"].idxmin()

    model_df = df[model_min_idx:model_min_idx + 1]
    loop_df = df[loop_min_idx:loop_min_idx + 1]

    with open(LENGTH_PATH + pdb_name + LENGTH_FILE, 'r') as file:
        length = [file.readline()]
        model_df["length"] = length
        loop_df["length"] = length

    return np.array(model_df), np.array(loop_df)


def extract_box_graphs():
    models, loops = generate_final_scores(sys.argv[1])

    generate_rmsd_graph(models, "Model")
    generate_rmsd_graph(loops, "Loop")


def extract_energy_vs_rmsd(pdb_path, pdb_name, score_type):
    df = pd.read_csv(pdb_path + "/" + pdb_name + "/scores.txt", sep=" ", header=None, names=COL)
    for i in range(5):
        loop_min_idx = df[df['type'] == "LOOP"][score_type].idxmin()
        model_min_idx = df[df['type'] == "MODEL"][score_type].idxmin()



if __name__ == '__main__':
    extract_box_graphs()


