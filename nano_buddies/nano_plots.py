import pandas as pd
import numpy as np
import os, sys
import enum
import matplotlib.pyplot as plt
from plotnine import *
import argparse


COL = ["type", "name", "dope_score", "soap_score", "rmsd"]
FULL_COL = ["type", "name", "dope_score", "soap_score", "rmsd", "length"]

LENGTH_PATH = "/cs/labs/dina/tomer.cohen13/lengths"
LENGTH_FILE = "nano_length.txt"

TOP_SCORES_N = 10


# Using enum class create enumerations
class ColInfo(enum.Enum):
    NAME = 1
    D_SCORE = 3
    S_SCORE = 5
    RMSD = 7
    LENGTH = 8


def get_scores_data(pdb_folder):
    """
    reads the scores.txt from the pdb_folder into df (column names = COL)
    :param pdb_folder: path of the pdb folder
    :return: df (len(COL) columns)
    """
    return pd.read_csv(os.path.join(pdb_folder, "scores.txt"), sep=" ", header=None, names=COL, usecols=[0,1,3,5,7])


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

    plot = ggplot(data,aes(x="factor(length)", y="rmsd")) + geom_boxplot(color="blue") + \
           geom_jitter(alpha=0.7, color="skyblue") + \
           ggtitle("rmsd VS. length for " + name + " results") + labs(x="length", y="RMSD")
    plot.save("rmsd_boxplot_" + name)


    # axes_boxplot = data.boxplot(column=["rmsd"], by=["length"], grid=False, showfliers=False)
    # axes_boxplot.set_title("rmsd VS. length for " + name + " results")
    # plt.savefig("rmsd_boxplot_" + name)


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
    df = get_scores_data(pdb_folder)

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


def get_min_rmsd_by_score(df, score_name):
    """
    returns the min rmsd index of the best 10 models according to score_name (lowest score), and the min rmsd index
    of all the models
    :param df: df to get the mins from
    :param score_name: score to use (dope_score/soap_score)
    :return: 2 floats
    """
    top_scores = pd.DataFrame.sort_values(df, by=score_name)[0:TOP_SCORES_N]
    top_min_index = top_scores["rmsd"].idxmin()

    rmsd_min_index = df["rmsd"].idxmin()

    return top_min_index, rmsd_min_index


def plot_points_one_pdb(pdb_folder, score_name):
    """
    plots a graph that shows the rmsd of the nanobody against its
    score
    :param pdb_path: pdb directory path
    :param score_name: score to use (dope_score/soap_score)
    :return: None
    """

    df = get_scores_data(pdb_folder)
    top_scores_index, rmsd_min_index = get_min_rmsd_by_score(df, score_name)

    top_scores_min = df.iloc[top_scores_index]["rmsd"]
    rmsd_min = df.iloc[rmsd_min_index]["rmsd"]

    # for legend names...
    if top_scores_index == rmsd_min_index:
        df.at[top_scores_index, "type"] = "MIN RMSD OF " + str(TOP_SCORES_N) + " BEST\nSCORES = MIN RMSD (" + \
                                          df.iloc[rmsd_min_index]["type"] + ")\nrmsd = %.2f" % rmsd_min
    else:
        df.at[top_scores_index, "type"] = "MIN RMSD OF " + str(TOP_SCORES_N) + " BEST\nSCORES (" + \
                                          df.iloc[top_scores_index]["type"] + ")\nrmsd = %.2f" % top_scores_min
        df.at[rmsd_min_index, "type"] = "MIN RMSD (" + df.iloc[rmsd_min_index]["type"] + ")\nrmsd = %.2f" % rmsd_min

    # create the plot
    pdb_name = os.path.basename(pdb_folder)
    y_label = score_name.replace("_", " ").upper()
    plot = ggplot(df) + geom_point(aes(x="rmsd", y=score_name, color="type"), alpha=0.7) \
           + ggtitle("RMSD agains " + y_label + " for nanobody " + pdb_name) + labs(x="RMSD", y=y_label) + \
           geom_text(aes(y=df.iloc[top_scores_index][score_name], x=top_scores_min, label="%.2f" % top_scores_min), size=8, ha="right", va="top") + \
           geom_text(aes(y=df.iloc[rmsd_min_index][score_name], x=rmsd_min, label="%.2f" % rmsd_min), size=8, ha="right", va="top")

    plot.save(pdb_name + "_" + score_name + "_plot")


def summery_rmsd_scores(directory, score):
    """
    saves a summery of all the min rmsd of each nanobody in the directory (both min rmsd and min rmsd of top 10 models
    according to score), save the data in the file directory/summery_rmsd_scores.csv
    :param directory: path to the pdbs folder
    :param score: score to use(dope_score,soap_score)
    :return: None
    """
    with open(os.path.join(directory, "summery_rmsd_scores.csv"), 'w') as output_file:
        first_pdb = True
        for file in os.listdir(directory):
            pdb_folder = os.path.join(directory, file)
            if os.path.isdir(pdb_folder) and os.path.exists(os.path.join(pdb_folder, "scores.txt")):
                one_pdb_rmsd_scores(pdb_folder, score).to_csv(output_file, header=first_pdb, index=False)
                first_pdb = False


def one_pdb_rmsd_scores(pdb_folder, score_name):
    """
    returns a df with the summery of one pdb nanobody (min rmsd by the best 10 models according to score_name, min rmsd
    of all models, typs, etc.)
    :param pdb_folder: the folder path of the pdb
    :param score_name: the score to use (dope_score, soap_score), string
    :return: df (1 row, 8 columns)
    """

    df = get_scores_data(pdb_folder)
    top_scores_index, rmsd_min_index = get_min_rmsd_by_score(df, score_name)
    top_scores_min = df.iloc[top_scores_index]["rmsd"]
    rmsd_min = df.iloc[rmsd_min_index]["rmsd"]

    pdb_name = os.path.basename(pdb_folder)
    output = pd.DataFrame({"PDB": [pdb_name], "SCORE_BEST_10": [df.iloc[top_scores_index][score_name]], "RMSD_BEST_10": [df.iloc[top_scores_index]["rmsd"]],
                           "TYPE_BEST_10": [df.iloc[top_scores_index]["type"]], "SCORE_MIN_RMSD": [df.iloc[rmsd_min_index][score_name]],
                           "MIN_RMSD": [df.iloc[rmsd_min_index]["rmsd"]], "TYPE_MIN_RMSD": [df.iloc[rmsd_min_index]["type"]],
                           "DIFF_RMSD": [top_scores_min - rmsd_min]})
    return output


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

    for file in os.listdir(directory):
        if n == 0:
            return
        pdb_folder = os.path.join(directory, file)
        if os.path.isdir(pdb_folder) and os.path.exists(os.path.join(pdb_folder, "scores.txt")):
            plot_points_one_pdb(pdb_folder, score)
            n -= 1


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="dirctory path containing the pdb files")
    parser.add_argument("score", help="score to use (dope_score/soap_score)")
    parser.add_argument("-b", "--boxplot", help="saves box plots", action="store_true")
    parser.add_argument("-p", "--points", help="saves points plots, gets number of pdbs to plot", type=int)

    args = parser.parse_args()
    if args.boxplot:  # if we want to create boxplot
        extract_box_graphs(args.directory)
    if args.points:  # if we want to create point plots (rmsd vs score)
        plot_rmsd_vs_score(args.directory, args.points, args.score)

    summery_rmsd_scores(args.directory, args.score)  # saves summery into csv file



