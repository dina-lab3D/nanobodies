import subprocess
import os
import numpy as np
import argparse
import re
import pandas as pd
from tqdm import tqdm
from plotnine import *


BAD_PDBS = ["4LAJ_1", "5J56_1", "5TOK_1", "5U4L_1"]
PLOTS_PATH = "/cs/labs/dina/tomer.cohen13/nanobodies/nano_buddies/docking_plots"
TOP_SCORES_N = 100


def make_data_one(pdb_name):

    #  get the rmsd column
    rmsd = subprocess.run("grep \"|\" docking_" + pdb_name + ".res | cut -d '|' -f2", shell=True, capture_output=True, universal_newlines=True).stdout

    #  get the soap_score column
    soap_score = subprocess.run("grep \"|\" soap_score_" + pdb_name + ".res | cut -d '|' -f2", shell=True, capture_output=True, universal_newlines=True).stdout

    #  get the transformation column
    trans = subprocess.run("grep \"|\" soap_score_" + pdb_name + ".res | cut -d '|' -f7", shell=True, capture_output=True, universal_newlines=True).stdout

    interface_rmsd = np.char.strip(np.char.strip(np.array(re.split(".+\(", rmsd))), ")")[1:] # split to array
    ligand_rmsd = np.char.strip(np.array(re.split(" rmsd   \n| \([0-9\.]+\)  \n", rmsd)))[1:-1]  # split to array
    soap_score = soap_score.replace(" ","").split("\n")[1:-1]  # split to array
    trans = trans.split("\n")[1:-1]
    indexes = list(range(1, len(soap_score) + 1))

    return ligand_rmsd.tolist(), interface_rmsd.tolist(), soap_score, trans, indexes


def make_data(folder):

    scores = []
    ligand = []
    interface = []
    names = []
    trans = []
    indexes = []
    for pdb_file in os.listdir(folder):
        #  loop/ model nanobody pdb
        if (pdb_file.startswith("model") or pdb_file.startswith("loop")) and pdb_file.endswith(".pdb") and "tr" not in pdb_file:  # TODO - change the if condition (less ugly...)
            pdb_name = pdb_file.split(".")[0]
            results = make_data_one(pdb_name)

            ligand += results[0]
            interface += results[1]
            scores += results[2]
            trans += results[3]
            indexes += results[4]
            names += len(results[0]) * [pdb_name]

    pd.DataFrame({"index": indexes, "names": names, "ligand_rmsd": ligand, "interface_rmsd": interface, "soap_score": scores, "trans": trans}).to_csv("dock_data.csv", index=False, header=True)


def dock_analyze(folder, to_plot, to_make_data):
    """

    :param folder:
    :param to_plot:
    :return:
    """
    os.chdir(folder)

    if to_make_data or not os.path.exists("dock_data.csv"):
        make_data(os.getcwd())

    data_df = pd.read_csv("dock_data.csv")

    data_df["ligand_rmsd"] = np.asfarray(data_df["ligand_rmsd"], float)  # cast to float
    data_df["interface_rmsd"] = np.asfarray(data_df["interface_rmsd"], float)  # cast to float
    data_df["soap_score"] = np.asfarray(data_df["soap_score"], float)  # cast to float

    top_scores = pd.DataFrame.sort_values(data_df, by="soap_score")[0:TOP_SCORES_N]
    min_score_idx = top_scores["ligand_rmsd"].idxmin()
    min_rmsd_idx = data_df["ligand_rmsd"].idxmin()

    ref_df = pd.read_csv("ref_scores.csv", header=None)
    ref_score = list(ref_df[ref_df[0] == "ref"][1])
    loop_model_score = list(ref_df[ref_df[0] == data_df.iloc[min_rmsd_idx]["names"]][1])

    df = pd.DataFrame({"..": ["        "], "min_score_name": data_df.iloc[min_score_idx]["names"] + " (" + str(min_score_idx) + ")", "min_score_soap": data_df.iloc[min_score_idx]["soap_score"], "min_score_ligand_rmsd": data_df.iloc[min_score_idx]["ligand_rmsd"],
                       "min_score_interface_rmsd": data_df.iloc[min_score_idx]["interface_rmsd"], "...": ["        "], "min_rmsd_name": data_df.iloc[min_rmsd_idx]["names"] + " (" + str(min_rmsd_idx) + ")", "min_rmsd_soap": data_df.iloc[min_rmsd_idx]["soap_score"],
                       "min_rmsd_ligand_rmsd": data_df.iloc[min_rmsd_idx]["ligand_rmsd"], "min_rmsd_interface_rmsd": data_df.iloc[min_rmsd_idx]["interface_rmsd"],
                       "....": ["        "], "diff_rmsd": [data_df.iloc[min_score_idx]["ligand_rmsd"] - data_df.iloc[min_rmsd_idx]["ligand_rmsd"]], "ref_score": ref_score,
                       "ref_model_loop_score": loop_model_score})

    if to_plot:
        plot_points(top_scores, os.path.basename(directory), data_df.iloc[min_rmsd_idx]["ligand_rmsd"], data_df.iloc[min_rmsd_idx]["soap_score"], min_score_idx, ref_score[0], loop_model_score[0])

    df.to_csv("dock_summery_best_rmsd_of_" + str(TOP_SCORES_N) + "_best_scores_" + os.path.basename(folder) + ".csv", header=True, index=False)
    os.chdir("..")
    return df


def plot_points(points_df, pdb_name, min_rmsd, min_rmsd_score, min_score_idx, ref_score, loop_model_score):

    ref_df = pd.read_csv("ref_scores.csv", header=None)
    plot = ggplot(points_df) + geom_point(aes(x="ligand_rmsd", y="soap_score", color="names"), alpha=0.7) + ggtitle("RMSD against Soap Score for " + str(TOP_SCORES_N) + " best dock results") + \
           geom_vline(xintercept=min_rmsd, linetype='dotted', color="red") + geom_hline(yintercept=min_rmsd_score, linetype='dotted', color="red")  + \
           geom_text(aes(y=points_df.loc[min_score_idx]["soap_score"], x=points_df.loc[min_score_idx]["ligand_rmsd"], label="%.2f" % points_df.loc[min_score_idx]["ligand_rmsd"]), size=8, ha="right", va="top") + \
           geom_text(aes(y=min_rmsd_score, x=min_rmsd, label="%.2f" % min_rmsd), size=8, ha="right", va="top") + geom_hline(yintercept=ref_score, linetype='dotted', color="green") + \
           geom_hline(yintercept=loop_model_score, linetype='dotted', color="blue") + \
           geom_text(aes(y=ref_score, x=1), size=8, ha="left", va="top", label="ref score") + \
           geom_text(aes(y=loop_model_score, x=1), size=8, ha="left", va="top", label="aligned loop or model score")

    plot.save(os.path.join(PLOTS_PATH, str(TOP_SCORES_N) + " rmsd_vs_soap_" + pdb_name))


def plots():

    dock_df = pd.read_csv("dock_summery_best_rmsd_of_" + str(TOP_SCORES_N) + "_best_scores" + ".csv")
    plot1 = ggplot(dock_df) + geom_point(aes(x="min_score_ligand_rmsd" , y="min_score_soap")) + \
            ggtitle("ligand RMSD of the dock result with minimum soap score against soap score") + labs(x="RMSD", y="Soap Score")
    plot2 = ggplot(dock_df) + geom_point(aes(x="min_rmsd_ligand_rmsd" , y="min_rmsd_soap")) + \
            ggtitle("ligand RMSD of the dock result with minimum RMSD against soap score") + labs(x="RMSD", y="Soap Score")

    plot1.save(os.path.join(PLOTS_PATH, str(TOP_SCORES_N) + "_min_score_plot"))
    plot2.save(os.path.join(PLOTS_PATH, "min_rmsd_plot"))


if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="directory path containing the pdb directories")
    parser.add_argument("-p", "--plots", help="plots the summery data (dock_data.csv has to exist first!)", action="store_true")
    parser.add_argument("-d", "--data", help="extract the data to dock_data.csv file", action="store_true")
    args = parser.parse_args()
    os.chdir(args.directory)
    first_pdb = True
    with open("dock_summery_best_rmsd_of_" + str(TOP_SCORES_N) + "_best_scores" + ".csv", 'w') as sum_file:
        for directory in tqdm(os.listdir(args.directory)):
            #  if the folder is pdb folder
            if os.path.isdir(directory) and re.fullmatch("[a-zA-Z0-9]{4}_[0-9]", directory):
                if directory in BAD_PDBS:  # TODO - delete when they are fixed
                    continue
                pdb_df = dock_analyze(os.path.join(args.directory, directory), args.plots, args.data)
                pdb_df.insert(0, column="pdb_name", value=[directory])
                pdb_df.to_csv(sum_file, header=first_pdb, index=False)
                first_pdb = False

    if args.plots:
        plots()











