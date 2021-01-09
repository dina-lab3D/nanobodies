import argparse
import os
import numpy as np
import pandas as pd
import pickle
from tqdm import tqdm
import re
from NanoNetUtils import generate_label, generate_input, valid_pdb


if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="dirctory path containing the pdb files")
    args = parser.parse_args()
    os.chdir(args.directory)
    failed_pdbs = pd.DataFrame(columns=["PDB", "FOLDER"])
    feature_matrix = []
    pdb_names = []
    input_matrix = []
    for directory in os.listdir(os.getcwd()):
        if os.path.isdir(directory) and re.fullmatch('[0-9]+', directory):  # directories 1,2,3,4...
            os.chdir(directory)
            print(directory)
            for pdb in tqdm(os.listdir(os.getcwd())):
                if os.path.isdir(pdb):
                    if not valid_pdb(pdb):
                        print(directory + ": " + pdb + ", FAILED")
                        failed_pdbs = failed_pdbs.append(pd.DataFrame({"PDB":[pdb], "FOLDER":[directory]}))
                        continue
                    feature_matrix.append(generate_label(pdb))
                    input_matrix.append(generate_input(pdb))
                    pdb_names.append(os.path.join(directory,pdb))
            os.chdir("..")
    feature_matrix = np.stack(feature_matrix, axis=0)
    input_matrix = np.stack(input_matrix, axis=0)
    labels_file_name = "nn_labels"
    input_file_name = "nn_input"
    pdb_names_file = "pdb_names"
    # if TEST:
    #     labels_file_name += "_test"
    #     pdb_names_file += "_test"
    #     input_file_name += "_test"
    # if BINS:
    #     labels_file_name += "_bins"
    #     pdb_names_file += "pdb_names"

    pickle.dump(feature_matrix, open(labels_file_name + ".pkl", "wb"))
    pickle.dump(input_matrix, open(input_file_name + ".pkl", "wb"))
    pickle.dump(np.array(pdb_names), open(pdb_names_file + ".pkl", "wb"))
    failed_pdbs.to_csv("nn_failed_pdbs.csv")