import argparse
import os
import subprocess
import re
import pandas as pd

INTERFACE = "/cs/labs/dina/tomer.cohen13/nanobodies/epiDock/interface "
HEADER = "/cs/usr/tomer.cohen13/lab/nanobodies/COVID_19/antibodies/anti_header.csv"


# def get_interface(folder):
#     os.chdir(folder)
#     print(folder)
#     df = None
#     if os.path.exists("nanobody_trans_0.pdb"):
#         subprocess.run(INTERFACE + " A H 6 nanobody_trans_*.pdb", shell=True)
#         df = pd.read_csv("epi.csv", header=None)
#     os.chdir("..")
#     return df
#
#
# if __name__ == '__main__':
#     parser = argparse.ArgumentParser()
#     parser.add_argument("directory", help="directory path containing the pdb directories")
#     args = parser.parse_args()
#     os.chdir(args.directory)
#
#     with open("epiDock_summey.csv", "w") as epi_file:
#         header_df = pd.read_csv(HEADER)
#         header_df.insert(0, column="PDB", value=None)
#         header_df.to_csv(epi_file, header=True, index=False)
#         for directory in os.listdir(os.getcwd()):
#             #  if the folder is pdb folder
#             if os.path.isdir(directory) and re.fullmatch("[a-zA-Z0-9]{4}_[0-9]", directory):
#                 epi_df = get_interface(directory)
#                 if epi_df is not None:
#                     epi_df.insert(0, column="PDB", value=directory)
#                     epi_df.to_csv(epi_file, header=False, index=False)


def get_interface(pdb, pdbs_df):
    dir_name = pdb.split("pdb")[1].split(".")[0]
    os.mkdir(dir_name)
    subprocess.run("mv " + pdb + " " + os.path.join(dir_name, pdb), shell=True)
    os.chdir(dir_name)

    spike_chains = pdbs_df[pdbs_df["PDB ID"] == dir_name.upper()]["SPIKE"].iloc[0]
    antibody_chains = pdbs_df[pdbs_df["PDB ID"] == dir_name.upper()]["ANTIBODY"].iloc[0].split(" ")
    antibody_names = pdbs_df[pdbs_df["PDB ID"] == dir_name.upper()]["NAME"].iloc[0].split(" ")

    dfs = []
    names = []
    for i in range(len(antibody_chains)):
        subprocess.run(INTERFACE + spike_chains + " " + antibody_chains[i] + " 6 " + pdb, shell=True)
        if os.path.exists("epi.csv"):
            anti_df = pd.read_csv("epi.csv", header=None).drop(1000, axis=1).astype(int).clip(0,1)
            dfs.append(anti_df)
            names.append(antibody_names[i])
            subprocess.run("mv epi.csv epi_" + antibody_names[i] + ".csv", shell=True)

    os.chdir("..")
    return dfs, names


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="directory path containing the pdb directories")
    parser.add_argument("csv", help="csv file with chains")
    args = parser.parse_args()
    os.chdir(args.directory)
    chains_df = pd.read_csv(args.csv)

    with open("epiDock_summey.csv", "w") as epi_file:
        header_df = pd.read_csv(HEADER)
        header_df.insert(0, column="PDB", value=None)
        header_df.insert(0, column="ANTIBODY", value=None)
        header_df.to_csv(epi_file, header=True, index=False)
        for pdb in os.listdir(os.getcwd()):
            #  if the folder is pdb folder
            if pdb.endswith(".ent"):
                epi_df, antibody_names = get_interface(pdb, chains_df)
                for df, name in zip(epi_df, antibody_names):
                    df.insert(0, column="ANTIBODY", value=name)
                    df.insert(0, column="PDB", value=pdb)
                    df.to_csv(epi_file, header=False, index=False)
