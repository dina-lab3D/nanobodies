
import argparse
import os
import re
import pandas as pd
import numpy as np
from tqdm import tqdm

MAX_N_TEMPLATES = 10


def get_seq(folder, rmsd_df):

    os.chdir(folder)
    pdb_name = os.path.basename(folder)

    seq_temp = np.array(pd.read_csv("seq_identity_filter",squeeze=True, header=None, skipinitialspace=True))
    df = pd.DataFrame(columns=["seq0", "seq1","seq2","seq3","seq4","seq5","seq6","seq7","seq8","seq9"])
    row = np.concatenate((seq_temp, (MAX_N_TEMPLATES-len(seq_temp)) * ["-"]))

    df.loc[0] = row
    df.insert(0, column="pdb_name", value=[pdb_name])
    df.insert(11, column="mean", value= np.mean(seq_temp))
    df.insert(12, column="median", value= np.median(seq_temp))
    df.insert(13, column="std", value= np.std(seq_temp))
    df.insert(14, column="rmsd", value=rmsd_df.iloc[0, 1])
    os.chdir("..")

    return df


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="directory path containing the pdb directories")
    args = parser.parse_args()
    os.chdir(args.directory)

    first_pdb = True
    output_file = "seq_indentity_summery.csv"
    rmsd_df = pd.read_csv(os.path.join(os.getcwd(), "model_summery", "summery_m5_l5_rmsd_soap_score.csv"), usecols=[0,2])
    # dock summery
    with open(output_file , 'w') as sum_file:
        for directory in tqdm(os.listdir(os.getcwd())):
            #  if the folder is pdb folder
            if os.path.isdir(directory) and re.fullmatch("[a-zA-Z0-9]{4}_[0-9]", directory):
                pdb_rmsd = rmsd_df[rmsd_df["PDB"] == directory]
                if pdb_rmsd.empty:
                    continue
                indentity_df = get_seq(directory, rmsd_df[rmsd_df["PDB"] == directory])
                indentity_df.to_csv(sum_file, header=first_pdb, index=False)
                first_pdb = False
