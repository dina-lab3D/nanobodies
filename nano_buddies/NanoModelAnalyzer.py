import pandas as pd
import numpy as np
import os
import enum
from plotnine import *
import argparse

# columns names for the data frame (of scores.txt)
COL = ["type", "name", "dope_score", "soap_score", "rmsd", "cdr1_rmsd", "cdr2_rmsd", "cdr3_rmsd"]
FULL_COL = ["type", "name", "dope_score", "soap_score", "rmsd", "cdr1_rmsd", "cdr2_rmsd", "cdr3_rmsd", "length"]

# columns names for file cdrs_dist
CDRS_COL = ["-3", "-2", "-1", "0 start", "1", "2", "3", "0 end"]

# path to the directory containing all the lengths of the pdbs
LENGTH_PATH = "/cs/labs/dina/tomer.cohen13/lengths"

# length file name
LENGTH_FILE = "nano_length.txt"


PLOTS_PATH = "model_plots"

# number of models to take in general
TOP_SCORES_N = 10

# number of models to take from loops
TOP_LOOP_N = 5

# number of models to take from models
TOP_MODEL_N = 5


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
    return pd.read_csv(os.path.join(pdb_folder, "scores.txt"), sep=" ", header=None, names=COL, usecols=[0,1,3,5,7,9,11,13])


def generate_rmsd_graph(folder, data, name):
    """
    saves the box plots of the length against minimal rmsd for all the
    nanobodies in data (box plot for each length)
    :param folder: pdbs folder path
    :param data: DataFrame with column names  = FULL_COL
    :param name: "Model"/"Loop"
    :return: None
    """
    data["rmsd"] = pd.to_numeric(data["rmsd"])  # cast to float
    data["length"] = pd.to_numeric(data["length"])  # cast to float

    plot = ggplot(data, aes(x="factor(length)", y="rmsd")) + geom_boxplot(color="black") + \
           geom_jitter(alpha=0.7, color="black", size=0.2) + \
           ggtitle("RMSD VS. length for " + name + " results") + labs(x="length", y="RMSD")
    plot.save(os.path.join(folder, PLOTS_PATH, "rmsd_boxplot_" + name), dpi=1000)


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

    generate_rmsd_graph(folder_path, models, "Model")
    generate_rmsd_graph(folder_path, loops, "Loop")


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
    cdr3_min_index = df["cdr3_rmsd"].idxmin()

    return top_min_index, rmsd_min_index, cdr3_min_index


def get_min_rmsd_by_type_score(df, score_name):
    """
    returns the min rmsd index of the best 10 models according to score_name (lowest score), and the min rmsd index
    of all the models
    :param df: df to get the mins from
    :param score_name: score to use (dope_score/soap_score)
    :return: 2 floats
    """
    top_loop_scores = pd.DataFrame.sort_values(df[df["type"] == "LOOP"], by=score_name)[0:TOP_LOOP_N]
    top_model_scores = pd.DataFrame.sort_values(df[df["type"] == "MODEL"], by=score_name)[0:TOP_MODEL_N]

    min_index = pd.concat([top_loop_scores, top_model_scores])["rmsd"].idxmin()
    rmsd_min_index = df["rmsd"].idxmin()
    cdr3_min_index = df["cdr3_rmsd"].idxmin()

    return min_index, rmsd_min_index, cdr3_min_index


def plot_points_one_pdb(folder, pdb_folder, score_name):
    """
    plots a graph that shows the rmsd of the nanobody against its
    score
    :param folder: pdbs folder
    :param pdb_folder: folder of the pdb
    :param score_name: score to use (dope_score/soap_score)
    :return: None
    """

    df = get_scores_data(pdb_folder)
    top_scores_index, rmsd_min_index, cdr3_min_index = get_min_rmsd_by_score(df, score_name)

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

    plot.save(os.path.join(folder, PLOTS_PATH, pdb_name + "_" + score_name + "_plot"), dpi=1000)


def plot_summery_scores(folder, input_file, score):

    df = pd.read_csv(input_file)
    name = os.path.basename(input_file).split(".")[0]

    left_of_line_best_10_rmsd= len(df[df['RMSD_BEST_10'] < 2]["RMSD_BEST_10"])
    left_of_line_min_rmsd = len(df[df['MIN_RMSD'] < 2]["MIN_RMSD"])

    # rmsd
    plot = ggplot(df) + geom_point(aes(x="RMSD_BEST_10", y="SCORE_BEST_10", color="TYPE_BEST_10"), alpha=0.7) + \
           geom_vline(xintercept=2, linetype='dotted', color="red", show_legend=True) + labs(x="RMSD", y="Soap score")+\
           ggtitle("RMSD of best score vs " + score.replace("_", " ") + " (left of line=" + str(left_of_line_best_10_rmsd) + ")")

    plot2 = ggplot(df) + geom_point(aes(x="MIN_RMSD", y="SCORE_MIN_RMSD", color="TYPE_MIN_RMSD"), alpha=0.7) + \
           geom_vline(xintercept=2, linetype='dotted', color="red", show_legend=True) + labs(x="RMSD", y="Soap Score")+\
            ggtitle("Min RMSD vs " + score.replace("_", " ") + " (left of line=" + str(left_of_line_min_rmsd) + ")")

    plot.save(os.path.join(folder, PLOTS_PATH, name), dpi=1000)
    plot2.save(os.path.join(folder, PLOTS_PATH, "summery_min_rmsd_" + score), dpi=1000)

    # cdrs
    for i in ["1","2","3"]:

        left_of_line_best_10_rmsd_cdr= len(df[df['CDR' + i] < 2]["CDR" + i])
        left_of_line_min_rmsd_cdr = len(df[df['MIN_CDR' + i] < 2]["MIN_CDR" + i])

        plot3 = ggplot(df) + geom_point(aes(x="CDR" + i, y="SCORE_BEST_10", color="TYPE_BEST_10"), alpha=0.7) + \
               geom_vline(xintercept=2, linetype='dotted', color="red", show_legend=True) + labs(x="RMSD", y="Soap score")+ \
               ggtitle("CDR" + i + " RMSD of best score vs " + score.replace("_", " ") + " (left of line=" + str(left_of_line_best_10_rmsd_cdr) + ")")

        plot4 = ggplot(df) + geom_point(aes(x="MIN_CDR" + i, y="SCORE_MIN_RMSD", color="TYPE_MIN_RMSD"), alpha=0.7) + \
            geom_vline(xintercept=2, linetype='dotted', color="red", show_legend=True) + labs(x="RMSD", y="Soap Score")+ \
            ggtitle("Min CDR" + i + " RMSD vs " + score.replace("_", " ") + " (left of line=" + str(left_of_line_min_rmsd_cdr) + ")")
        plot3.save(os.path.join(folder, PLOTS_PATH, name + '_cdr' + i), dpi=1000)
        plot4.save(os.path.join(folder, PLOTS_PATH, "summery_min_rmsd_" + score + "_cdr" + i), dpi=1000)


def summery_rmsd_scores(directory, score):
    """
    saves a summery of all the min rmsd of each nanobody in the directory (both min rmsd and min rmsd of top 10 models
    according to score), save the data in the file directory/summery_rmsd_scores.csv
    :param directory: path to the pdbs folder
    :param score: score to use(dope_score,soap_score)
    :return: None
    """
    summery_best_10 = os.path.join(directory, "model_summery", "summery_10_rmsd_" + score + ".csv")
    summery_best_5_by_type = os.path.join(directory, "model_summery", "summery_" + "m" + str(TOP_MODEL_N) + "_l" + str(TOP_LOOP_N) + "_rmsd_" + score + ".csv")

    with open(summery_best_10, 'w') as output_file_10:
        with open(summery_best_5_by_type, "w") as output_file_5_5:
            first_pdb = True
            for file in os.listdir(directory):
                pdb_folder = os.path.join(directory, file)
                if os.path.isdir(pdb_folder) and os.path.exists(os.path.join(pdb_folder, "scores.txt")):
                    one_pdb_rmsd_scores(pdb_folder, score, get_min_rmsd_by_score).to_csv(output_file_10, header=first_pdb, index=False)
                    one_pdb_rmsd_scores(pdb_folder, score, get_min_rmsd_by_type_score).to_csv(output_file_5_5, header=first_pdb, index=False)
                    first_pdb = False

    plot_summery_scores(directory, summery_best_10, score)
    plot_summery_scores(directory, summery_best_5_by_type, score)

    plot_boxplot_cdrs_rmsd(directory, summery_best_10)
    plot_boxplot_cdrs_rmsd(directory, summery_best_5_by_type)

    summery_differences(directory, score, summery_best_10, summery_best_5_by_type)


def plot_boxplot_cdrs_rmsd(folder, input_file):
    """

    :param folder:
    :param input_file:
    :param score:
    :return:
    """
    df = pd.read_csv(input_file, usecols=[2,3,4,5]).rename(columns={"RMSD_BEST_10": "ALL"}).melt(value_vars=["ALL","CDR1","CDR2","CDR3"])
    print(df)
    name = os.path.basename(input_file).split(".")[0]

    plot = ggplot(df, aes(x="factor(variable)", y="value")) + geom_boxplot(color="black") + \
           geom_jitter(alpha=0.7, color="black", size=0.2) + \
           ggtitle("RMSD of Nanobodies CDRs") + labs(x="Region", y="RMSD")
    plot.save(os.path.join(folder, PLOTS_PATH,name + "_cdrs_rmsd_boxplot"), dpi=1000)


def summery_differences(directory, score_name, summery_best_10_path, summery_best_5_by_type_path):
    """

    :param summery_best_10_path:
    :param summery_best_5_by_type_path:
    :param directory
    :param score_name
    :return:
    """
    df_10 = pd.read_csv(summery_best_10_path)
    df_5_5 = pd.read_csv(summery_best_5_by_type_path)

    rmsd_sum1, rmsd_sum2 = df_10["RMSD_BEST_10"].sum(), df_5_5["RMSD_BEST_10"].sum()
    diff_sum1, diff_sum2 = df_10["DIFF_RMSD"].sum(), df_5_5["DIFF_RMSD"].sum()
    diff_03_1, diff_03_2 = len(df_10[df_10["DIFF_RMSD"] < 0.3]["DIFF_RMSD"]), len(df_5_5[df_5_5["DIFF_RMSD"] < 0.3]["DIFF_RMSD"])
    rmsd_under_2_1, rmsd_under_2_2 = len(df_10[df_10["RMSD_BEST_10"] < 2]["RMSD_BEST_10"]), len(df_5_5[df_5_5["RMSD_BEST_10"] < 2]["RMSD_BEST_10"])
    rmsd_under_15_1, rmsd_under_15_2 = len(df_10[df_10["RMSD_BEST_10"] < 1.5]["RMSD_BEST_10"]), len(df_5_5[df_5_5["RMSD_BEST_10"] < 1.5]["RMSD_BEST_10"])
    rmsd_under_1_1, rmsd_under_1_2 = len(df_10[df_10["RMSD_BEST_10"] < 1]["RMSD_BEST_10"]), len(df_5_5[df_5_5["RMSD_BEST_10"] < 1]["RMSD_BEST_10"])
    rmsd_under_05_1, rmsd_under_05_2 = len(df_10[df_10["RMSD_BEST_10"] < 0.5]["RMSD_BEST_10"]), len(df_5_5[df_5_5["RMSD_BEST_10"] < 0.5]["RMSD_BEST_10"])

    df1 = pd.DataFrame({"type": ["10"], "rmsd_sum": [rmsd_sum1], "diff_sum": [diff_sum1], "diff_num_under_0.3": [diff_03_1], "num_rmsd_under_2": [rmsd_under_2_1],
                       "num_rmsd_under_1.5": [rmsd_under_15_1], "num_rmsd_under_1": [rmsd_under_1_1], "num_rmsd_under_0.5": [rmsd_under_05_1]})
    df2 = pd.DataFrame({"type": ["5"], "rmsd_sum": [rmsd_sum2], "diff_sum": [diff_sum2], "diff_num_under_0.3": [diff_03_2], "num_rmsd_under_2": [rmsd_under_2_2],
                       "num_rmsd_under_1.5": [rmsd_under_15_2], "num_rmsd_under_1": [rmsd_under_1_2], "num_rmsd_under_0.5": [rmsd_under_05_2]})

    pd.concat([df1, df2]).to_csv(os.path.join(directory, "model_summery", "diff_summery_" + score_name + "m" + str(TOP_MODEL_N) + "l" + str(TOP_LOOP_N) + ".csv"), header=True, index=False)


def one_pdb_rmsd_scores(pdb_folder, score_name, min_func):
    """
    returns a df with the summery of one pdb nanobody (min rmsd by the best 10 models according to score_name, min rmsd
    of all models, typs, etc.)
    :param pdb_folder: the folder path of the pdb
    :param score_name: the score to use (dope_score, soap_score), string
    :param min_func: min function to use
    :return: df (1 row, 8 columns)
    """

    df = get_scores_data(pdb_folder)
    top_scores_index, rmsd_min_index, cdr3_min_index = min_func(df, score_name)
    top_scores_min = df.iloc[top_scores_index]["rmsd"]
    rmsd_min = df.iloc[rmsd_min_index]["rmsd"]
    cdr3_rmsd_min = df.iloc[cdr3_min_index]["cdr3_rmsd"]

    pdb_name = os.path.basename(pdb_folder)
    output = pd.DataFrame({"PDB": [pdb_name], "SCORE_BEST_10": [df.iloc[top_scores_index][score_name]], "RMSD_BEST_10": [df.iloc[top_scores_index]["rmsd"]],
                           "CDR1": [df.iloc[top_scores_index]["cdr1_rmsd"]], "CDR2": [df.iloc[top_scores_index]["cdr2_rmsd"]], "CDR3": [df.iloc[top_scores_index]["cdr3_rmsd"]],
                           "TYPE_BEST_10": [df.iloc[top_scores_index]["type"]], "SCORE_MIN_RMSD": [df.iloc[rmsd_min_index][score_name]],
                           "MIN_RMSD": [df.iloc[rmsd_min_index]["rmsd"]], "MIN_CDR1": [df.iloc[rmsd_min_index]["cdr1_rmsd"]], "MIN_CDR2": [df.iloc[rmsd_min_index]["cdr2_rmsd"]],
                           "MIN_CDR3": [df.iloc[rmsd_min_index]["cdr3_rmsd"]], "TYPE_MIN_RMSD": [df.iloc[rmsd_min_index]["type"]],
                           "DIFF_RMSD": [top_scores_min - rmsd_min], "MIN_CDR3_RMSD": cdr3_rmsd_min})
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
            plot_points_one_pdb(directory, pdb_folder, score)
            n -= 1


def summery_cdr3(directory):

    mean_file_name = os.path.join(directory, "model_summery", "cdr3_mean.csv")
    all_file_name = os.path.join(directory, "model_summery", "cdr3_all.csv")
    with open(mean_file_name, 'w') as mean_file:
        with open(all_file_name, "w") as all_file:
            first_pdb = True
            for file in os.listdir(directory):
                pdb_folder = os.path.join(directory, file)
                if os.path.isdir(pdb_folder) and os.path.exists(os.path.join(pdb_folder, "cdrs_dist")):
                    df = pd.read_csv(os.path.join(pdb_folder, "cdrs_dist"), sep=" ", header=None, names=CDRS_COL, usecols=[3,5,7,9,11,13,15,17])
                    df.insert(0, column="pdb", value=[os.path.basename(pdb_folder)]*len(df["-1"]))

                    df.to_csv(all_file, header=first_pdb, index=False)
                    mean_series = pd.concat([pd.Series(os.path.basename(pdb_folder)), df.mean(axis=0)])
                    pd.DataFrame(columns=mean_series.index, data=[mean_series]).to_csv(mean_file, header=first_pdb, index=False)
                    first_pdb = False

    plot = ggplot(pd.melt(pd.read_csv(mean_file_name), id_vars=['0'], value_vars=CDRS_COL), aes(x="factor(variable)", y="value")) + geom_boxplot(color="black") + \
           geom_jitter(alpha=0.3, color="black", size=0.2) + ggtitle("CDR3 frame vs distance from ref") + labs(x="cdr3 frame", y="distance from ref (Angstram)")
    plot.save(os.path.join(args.directory, PLOTS_PATH, "cdr3_frames_boxplot"), dpi=1000)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="directory path containing the pdb directories")
    parser.add_argument("score", help="score to use (dope_score/soap_score)")
    parser.add_argument("-b", "--boxplot", help="saves box plots", action="store_true")
    parser.add_argument("-p", "--points", help="saves points plots, gets number of pdbs to plot", type=int)
    parser.add_argument("-s", "--summery", help="saves sumery plots and creats summry csv", action="store_true")

    args = parser.parse_args()
    if not os.path.isdir(os.path.join(args.directory, PLOTS_PATH)):
        os.mkdir(os.path.join(args.directory, PLOTS_PATH))

    if args.boxplot:  # if we want to create boxplot
        extract_box_graphs(args.directory)
    if args.points:  # if we want to create point plots (rmsd vs score)
        plot_rmsd_vs_score(args.directory, args.points, args.score)
    if args.summery:  # if we want to create point plots (rmsd vs score)
        if not os.path.isdir(os.path.join(args.directory, "model_summery")):
            os.mkdir(os.path.join(args.directory, "model_summery"))
        summery_rmsd_scores(args.directory, args.score)  # saves summery into csv file
        summery_cdr3(args.directory)


