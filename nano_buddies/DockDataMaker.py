import os
import argparse
import subprocess
import pandas as pd


def make_data_one(pdb_name, xl):
    """
    organises the docking data of pdb_name
    :param pdb_name: the pdb to use
    :param xl: true if used cross links in the docking
    :return: DataFrame
    """
    rows_to_skip = 25  # intro length
    if xl:
        rows_to_skip = 30
    docking_df = pd.read_csv("docking_" + pdb_name + ".res", sep="|", header=0, usecols=[1], skipinitialspace=True, skiprows = list(range(0,rows_to_skip)), skipfooter=1)
    score_df = pd.read_csv("soap_score_" + pdb_name + ".res" , sep="|", header=0, usecols=[0, 1, 6], skipinitialspace=True, skiprows=list(range(0,3)))
    rmsds_df = docking_df.iloc[:, 0].str.extract("([\d\.]+) \(([\d\.]+)\)")

    return pd.DataFrame({"index": score_df.iloc[:, 0], "names": [pdb_name] * len(rmsds_df), "ligand_rmsd": rmsds_df.iloc[:, 0], "interface_rmsd": rmsds_df.iloc[:, 1], "soap_score": score_df.iloc[:, 1], "trans": score_df.iloc[:, 2]})


def combine_refs(folder):
    """
    organises all the ref docking data into a file ref_scores.csv
    :param folder: the folder containing the docking data
    :return:  None
    """
    names = []
    scores = []

    for file in os.listdir(folder):
        if file.startswith("no_trans") or file == "soap_score_ref_temp.res":
            soap_score = subprocess.run("grep \"|\" " + file + " | cut -d '|' -f2", shell=True, capture_output=True, universal_newlines=True).stdout
            soap_score = soap_score.replace(" ","").split("\n")[1:-1]  # split to array
            names.append(file.split(".")[0].split("score_")[-1])
            scores.append(soap_score[0])
            os.remove(file)
    df = pd.DataFrame({"name": names, "soap_score": scores})
    df.to_csv("ref_scores.csv", index=False, header=True)


def make_data(folder, xl):
    """
    organises all the docking data into a file named dock_data.csv
    :param folder: folder containing the pdb docking data
    :param xl: true if used cross links in the docking
    :return: None
    """
    first = True
    with open("dock_data.csv", 'w') as file:
        for pdb_file in os.listdir(folder):
            #  loop/ model nanobody pdb
            if pdb_file.startswith("params") and pdb_file != "params_ref.txt":  # goes over all the models used for docking
                pdb_name = (pdb_file.split("_")[1] + "_" + pdb_file.split("_")[2]).split(".")[0]
                df = make_data_one(pdb_name, xl)
                df.to_csv(file, header=first, index=False)
                first = False


if __name__ == '__main__':
    """
    organises all the docking data into a file named dock_data.csv and a file named ref_scores.csv
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="directory path to the pdb folder")
    parser.add_argument("-xl", "--crosslinks", help="if used cross links in the docking", action="store_true")
    args = parser.parse_args()
    os.chdir(args.directory)

    # the docking data of the ref pdbs
    combine_refs(os.getcwd())

    # all the docking data without the ref pdbs
    make_data(os.getcwd(), args.crosslinks)

    os.chdir("..")




