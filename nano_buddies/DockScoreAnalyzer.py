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


def ref_analyze(folder):

    os.chdir(folder)

    docking_df = pd.read_csv("pydock_results/docking_ref.res", sep="|", header=0, usecols=[1], skipinitialspace=True, skiprows = list(range(0,25)), skipfooter=1)
    score_df = pd.read_csv("pydock_results/soap_score_ref.res" , sep="|", header=0, usecols=[0, 1, 6], skipinitialspace=True, skiprows=list(range(0,3)))
    if (score_df.empty) or (docking_df.empty):
        os.chdir("..")
        return None
    rmsds_df = docking_df.iloc[:, 0].str.extract("([\d\.]+) \(([\d\.]+)\)")
    ref_df = pd.DataFrame({"names": [os.path.basename(folder)] * len(rmsds_df), "ligand_rmsd": rmsds_df.iloc[:, 0], "interface_rmsd": rmsds_df.iloc[:, 1], "soap_score": score_df.iloc[:, 1]})
    ref_df["ligand_rmsd"] = np.asfarray(ref_df["ligand_rmsd"], float)

    min_rmsd_idx = ref_df["ligand_rmsd"].idxmin()
    os.chdir("..")
    return ref_df[min_rmsd_idx:min_rmsd_idx + 1]


def dock_analyze(folder, to_plot, use_cluster, cross_links):
    """

    :param folder:
    :param to_plot:
    :return:
    """
    in_folder = "pydock_results"
    if cross_links:
        in_folder += "_xl_" + cross_links

    os.chdir(os.path.join(folder, in_folder))
    if use_cluster:
        data_df = pd.read_csv("soap_score_cluster.res", sep="|", header=None, skipinitialspace=True, skiprows=list(range(0, 2)), names=["number", "index", "names", "ligand_rmsd", "interface_rmsd", "soap_score", "cluster_size", "transformation"])
    else:
        try:
            data_df = pd.read_csv("dock_data.csv")
        except pd.errors.EmptyDataError:
            data_df = pd.DataFrame()
    if data_df.empty:
        os.chdir("..")
        os.chdir("..")
        return pd.DataFrame({"..": ["        "], "min_score_name": "-", "min_score_soap": "-", "min_score_ligand_rmsd": "-",
                             "min_score_interface_rmsd": "-", "...": ["        "], "min_rmsd_name": "-", "min_rmsd_soap": "-",
                             "min_rmsd_ligand_rmsd": "-", "min_rmsd_interface_rmsd": "-", "....": ["        "], "diff_rmsd": "-", "ref_score": "-",
                             "ref_model_loop_score": "-"})
    data_df["ligand_rmsd"] = np.asfarray(data_df["ligand_rmsd"], float)  # cast to float
    data_df["interface_rmsd"] = np.asfarray(data_df["interface_rmsd"], float)  # cast to float
    data_df["soap_score"] = np.asfarray(data_df["soap_score"], float)  # cast to float

    top_scores = pd.DataFrame.sort_values(data_df, by="soap_score")[0:TOP_SCORES_N]
    min_score_idx = top_scores["ligand_rmsd"].idxmin()
    min_rmsd_idx = data_df["ligand_rmsd"].idxmin()

    ref_df = pd.read_csv("ref_scores.csv", header=None)
    ref_score = list(ref_df[ref_df[0] == "ref"][1])
    loop_model_score = list(ref_df[ref_df[0] == data_df.iloc[min_rmsd_idx]["names"].strip()][1])

    df = pd.DataFrame({"..": ["        "], "min_score_name": data_df.iloc[min_score_idx]["names"] + " (" + str(min_score_idx+1) + ")", "min_score_soap": data_df.iloc[min_score_idx]["soap_score"], "min_score_ligand_rmsd": data_df.iloc[min_score_idx]["ligand_rmsd"],
                       "min_score_interface_rmsd": data_df.iloc[min_score_idx]["interface_rmsd"], "...": ["        "], "min_rmsd_name": data_df.iloc[min_rmsd_idx]["names"] + " (" + str(min_rmsd_idx+1) + ")", "min_rmsd_soap": data_df.iloc[min_rmsd_idx]["soap_score"],
                       "min_rmsd_ligand_rmsd": data_df.iloc[min_rmsd_idx]["ligand_rmsd"], "min_rmsd_interface_rmsd": data_df.iloc[min_rmsd_idx]["interface_rmsd"],
                       "....": ["        "], "diff_rmsd": [data_df.iloc[min_score_idx]["ligand_rmsd"] - data_df.iloc[min_rmsd_idx]["ligand_rmsd"]], "ref_score": ref_score,
                       "ref_model_loop_score": loop_model_score})    # TODO : change the index of the model with minimal score and the model with minimal rmsd (now its wrong)

    if to_plot:
        plot_points(top_scores, os.path.basename(directory), data_df.iloc[min_rmsd_idx]["ligand_rmsd"], data_df.iloc[min_rmsd_idx]["soap_score"], min_score_idx, float(ref_score[0]), float(loop_model_score[0]), use_cluster, cross_links)

    df.to_csv("dock_summery_best_rmsd_of_" + str(TOP_SCORES_N) + "_best_scores_" + os.path.basename(folder) + ".csv", header=True, index=False)
    os.chdir("..")
    os.chdir("..")
    return df


def plot_points(points_df, pdb_name, min_rmsd, min_rmsd_score, min_score_idx, ref_score, loop_model_score, use_cluster, cross_links):

    ref_df = pd.read_csv("ref_scores.csv", header=None)
    plot = ggplot(points_df) + geom_point(aes(x="ligand_rmsd", y="soap_score", color="names"), alpha=0.7) + ggtitle("RMSD against Soap Score for " + str(TOP_SCORES_N) + " best dock results") + \
           geom_vline(xintercept=min_rmsd, linetype='dotted', color="red") + geom_hline(yintercept=min_rmsd_score, linetype='dotted', color="red")  + \
           geom_text(aes(y=points_df.loc[min_score_idx]["soap_score"], x=points_df.loc[min_score_idx]["ligand_rmsd"], label="%.2f" % points_df.loc[min_score_idx]["ligand_rmsd"]), size=8, ha="right", va="top") + \
           geom_text(aes(y=min_rmsd_score, x=min_rmsd, label="%.2f" % min_rmsd), size=8, ha="right", va="top") + geom_hline(yintercept=ref_score, linetype='dotted', color="green") + \
           geom_hline(yintercept=loop_model_score, linetype='dotted', color="blue") + \
           geom_text(aes(y=ref_score, x=1), size=8, ha="left", va="top", label="ref score") + \
           geom_text(aes(y=loop_model_score, x=1), size=8, ha="left", va="top", label="aligned loop or model score")

    folder_name = os.path.join(PLOTS_PATH, "plots_" + str(TOP_SCORES_N))
    if use_cluster:
        folder_name += "_cluster"
    if cross_links:
        folder_name += "_xl_" + cross_links
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)
    plot_name = os.path.join(folder_name, str(TOP_SCORES_N) + "_rmsd_vs_soap_" + pdb_name)
    if use_cluster:
        plot_name += "_c"
    if cross_links:
        plot_name += "_xl_" + cross_links
    plot.save(plot_name, verbose=False)


def plots(use_cluster, cross_links):
    input_file = "dock_summery_best_rmsd_of_" + str(TOP_SCORES_N) + "_best_scores"
    if use_cluster:
        input_file += "_cluster"
    if cross_links:
        input_file += "_xl_" + cross_links
    dock_df = pd.read_csv(input_file + ".csv", na_values="-", na_filter=True)
    plot1 = ggplot(dock_df) + geom_point(aes(x="min_score_ligand_rmsd" , y="min_score_soap")) + \
            ggtitle("ligand RMSD of the dock result with minimum soap score against soap score") + labs(x="RMSD", y="Soap Score")
    plot2 = ggplot(dock_df) + geom_point(aes(x="min_rmsd_ligand_rmsd" , y="min_rmsd_soap")) + \
            ggtitle("ligand RMSD of the dock result with minimum RMSD against soap score") + labs(x="RMSD", y="Soap Score")
    plot1_name = os.path.join(PLOTS_PATH, str(TOP_SCORES_N) + "_min_score_plot")
    plot2_name = os.path.join(PLOTS_PATH, "min_rmsd_plot")
    if use_cluster:
        plot1_name += "_cluster"
        plot2_name += "_cluster"
    if cross_links:
        plot1_name += "_xl_" + cross_links
        plot2_name += "_xl_" + cross_links
    plot1.save(plot1_name, verbose=False)
    plot2.save(plot2_name, verbose=False)


if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="directory path containing the pdb directories")
    parser.add_argument("-p", "--plots", help="plots the summery data (dock_data.csv has to exist first!)", action="store_true")
    parser.add_argument("-c", "--cluster", help="use the data after clustering", action="store_true")
    parser.add_argument("-xl", "--crosslinks", help="use the data after cross links")
    args = parser.parse_args()
    os.chdir(args.directory)
    # first_pdb = True
    # output_file = "dock_summery_best_rmsd_of_" + str(TOP_SCORES_N) + "_best_scores"
    #
    # if args.cluster:
    #     output_file += "_cluster"
    # if args.crosslinks:
    #     output_file += "_xl_" + args.crosslinks
    #
    # with open(output_file + ".csv", 'w') as sum_file:
    #     for directory in tqdm(os.listdir(args.directory)):
    #         #  if the folder is pdb folder
    #         if os.path.isdir(directory) and re.fullmatch("[a-zA-Z0-9]{4}_[0-9]", directory):
    #             if directory in BAD_PDBS:  # TODO - delete when they are fixed
    #                 continue
    #             pdb_df = dock_analyze(os.path.join(args.directory, directory), args.plots, args.cluster, args.crosslinks)
    #             pdb_df.insert(0, column="pdb_name", value=[directory])
    #             pdb_df.to_csv(sum_file, header=first_pdb, index=False)
    #             first_pdb = False
    first_pdb = True
    with open("ref_docking_summery.csv", "w") as ref_sum_file:
        for directory in tqdm(os.listdir(args.directory)):
            #  if the folder is pdb folder
            if os.path.isdir(directory) and re.fullmatch("[a-zA-Z0-9]{4}_[0-9]", directory):
                if directory in BAD_PDBS:  # TODO - delete when they are fixed
                    continue

                ref_df = ref_analyze(os.path.join(args.directory, directory))
                if ref_df is None:
                    print(directory)
                    continue
                ref_df.to_csv(ref_sum_file, header=first_pdb, index=False)
                first_pdb = False
    if args.plots:
        plots(args.cluster, args.crosslinks)











