import pandas as pd
import numpy as np
import os, sys
import enum
import matplotlib.pyplot as plt
from plotnine import *
import argparse


COL = ["type", "name", "2", "dope_score", "4", "soap_score","6" , "rmsd"]
FULL_COL = ["type", "name", "2", "dope_score", "4", "soap_score","6" , "rmsd", "length"]

LENGTH_PATH = "/cs/labs/dina/tomer.cohen13/lengths"
LENGTH_FILE = "nano_length.txt"


# Using enum class create enumerations
class ColInfo(enum.Enum):
    NAME = 1
    D_SCORE = 3
    S_SCORE = 5
    RMSD = 7
    LENGTH = 8


def generate_rmsd_graph(data, name):
    """
    saves the box plots of the length against minimal rmsd for all the
    nanobodies in data (box plot for each length)
    :param data: DataFrame with column names  = FULL_COL
    :param name: "Model"/"Loop"
    :return: None
    """
    data["rmsd"] = pd.to_numeric(data["rmsd"])
    data["length"] = pd.to_numeric(data["length"])
    axes_boxplot = data.boxplot(column=["rmsd"], by=["length"], grid=False, showfliers=False)
    axes_boxplot.set_title("rmsd VS. length for " + name + " results")
    plt.savefig("rmsd_boxplot_" + name)



def generate_final_scores(folder_path):
    """
    returns 2 DataFrames (MODEL, LOOP) each with the minimal rmsd for each
    nanobody and its length
    :param folder_path: pdbs directory path
    :return: 2 DataFrames, size : (number of nanobodies in folder_path) * (len(FULL_COL))
    """
    models, loops= [], []
    # create a list of the scores from the pdb files, divided into base model and sampled loops
    for file in os.listdir(folder_path):
        pdb_folder = os.path.join(folder_path, file)
        if os.path.isdir(pdb_folder) and os.path.exists(os.path.join(pdb_folder, "scores.txt")):
            model, loop = parse_one_pdb(pdb_folder)
            models.append(model.tolist()[0])
            loops.append(loop.tolist()[0])
    return pd.DataFrame(data=models, columns=FULL_COL), pd.DataFrame(data=loops, columns=FULL_COL)


def parse_one_pdb(pdb_folder):
    """
    returns 2 numpy arrays (MODEL, LOOP) that contains the minimal rmsd of the
    nanobody in the pdb_folder and also the nanobody length (the indexes in the
    numpy array correspond to FULL_COL)
    :param pdb_folder:  pdb directory path
    :return: 2 numpy arrays, size: 1 * len(FULL_COL)
    """
    df = pd.read_csv(os.path.join(pdb_folder, "scores.txt"), sep=" ", header=None, names=COL)

    loop_min_idx = df[df['type'] == "LOOP"]["rmsd"].idxmin()
    model_min_idx = df[df['type'] == "MODEL"]["rmsd"].idxmin()

    model_df = df[model_min_idx:model_min_idx + 1]
    loop_df = df[loop_min_idx:loop_min_idx + 1]

    pdb_name = os.path.basename(pdb_folder)
    pdb_length_path = os.path.join(LENGTH_PATH, pdb_name, LENGTH_FILE)
    with open(pdb_length_path, 'r') as file:
        length = [file.readline()]
        model_df["length"] = length
        loop_df["length"] = length

    return np.array(model_df), np.array(loop_df)


def extract_box_graphs(folder_path):
    """
    plots the box graphs of the pdbs in the folder_path (x - length,
    y - minimal rmsd score)
    :param folder_path: pdbs directory path
    :return: None
    """
    models, loops = generate_final_scores(folder_path)

    generate_rmsd_graph(models, "Model")
    generate_rmsd_graph(loops, "Loop")


# def extract_energy_vs_rmsd(pdb_path, pdb_name, score_type):
#     """'
#     """
#     df = pd.read_csv(pdb_path + "/" + pdb_name + "/scores.txt", sep=" ", header=None, names=COL)
#     for i in range(5):
#         loop_min_idx = df[df['type'] == "LOOP"][score_type].idxmin()
#         model_min_idx = df[df['type'] == "MODEL"][score_type].idxmin()


def plot_points_one_pdb(pdb_folder, score_name):
    """
    plots a graph that shows the rmsd of the nanobody against its
    score
    :param pdb_path: pdb directory path
    :param score_name: score to use (dope_score/soap_score)
    :return: None
    """
    pdb_name = os.path.basename(pdb_folder)
    df = pd.read_csv(os.path.join(pdb_folder, "scores.txt"), sep=" ", header=None, names=COL)
    plot = ggplot(df) + geom_point(aes(x=score_name, y="rmsd"), color="blue") \
           + ggtitle("RMSD agains " + score_name + " for nanobody " + pdb_name)
    plot.save(pdb_name + "_" + score_name + "_plot")


def plot_rmsd_vs_score(directory, n, score):
    """
    plots n graphs for n different nanobodies that have a pdb directory in the
    "directory" argument, the graphs show the rmsd of the nanobody against its
    score
    :param directory: pdbs directory path
    :param n: number of nanobodies to plot (int)
    :param score: score to use (dope_score/soap_score)
    :return: None
    """
    while n != 0:
        for file in os.listdir(directory):
            pdb_folder = os.path.join(directory, file)
            if os.path.isdir(pdb_folder) and os.path.exists(os.path.join(pdb_folder, "scores.txt")):
                plot_points_one_pdb(pdb_folder, score)
                n -= 1


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="dirctory path containing the pdb files")
    parser.add_argument("n", type=int, help="number of plots to save (each for a different nanobody)")
    parser.add_argument("score", help="score to use (dope_score/soap_score)")
    parser.add_argument("-b", "--boxplot", help="saves box plots", action="store_true")

    args = parser.parse_args()
    if args.boxplot:
        extract_box_graphs(args.directory)

    plot_rmsd_vs_score(args.directory, args.n, args.score)  # rmsd vs score



