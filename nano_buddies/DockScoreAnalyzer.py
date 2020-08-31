import subprocess
import os
import numpy as np
import argparse
import re
import pandas as pd

BAD_PDBS = ["4LAJ_1", "5J56_1", "5TOK_1", "5U4L_1"]
def main(directory):
    """

    :param directory:
    :return:
    """
    os.chdir(directory)
    scores = []
    ligand = []
    interface = []
    names = []
    for pdb_file in os.listdir(os.getcwd()):
        #  loop/ model nanobody pdb
        if (pdb_file.startswith("model") or pdb_file.startswith("loop")) and pdb_file.endswith(".pdb") and "tr" not in pdb_file:  # TODO - change the if condition (less ugly...)
            pdb_name = pdb_file.split(".")[0]

            #  get the rmsd column
            rmsd = subprocess.run("grep \"|\" docking_" + pdb_name + ".res | cut -d '|' -f2", shell=True, capture_output=True, universal_newlines=True).stdout

            #  get the soap_score column
            soap_score = subprocess.run("grep \"|\" soap_score_" + pdb_name + ".res | cut -d '|' -f2", shell=True, capture_output=True, universal_newlines=True).stdout

            interface_rmsd = np.char.strip(np.char.strip(np.array(re.split(".+\(", rmsd))), ")")[1:-2]
            ligand_rmsd = np.char.strip(np.array(re.split(" rmsd   \n| \([0-9\.]+\)  \n", rmsd)))[1:-2]  # split to array
            soap_score = soap_score.replace(" ","").split("\n")[1:-2]  # split to array

            ligand += ligand_rmsd.tolist()
            interface += interface_rmsd.tolist()
            scores += soap_score
            names += len(soap_score) * [pdb_name]

    ligand = np.asfarray(ligand, float)  # cast to float
    interface = np.asfarray(interface, float)  # cast to float
    scores = np.asfarray(scores, float)  # cast to float

    min_score_idx = np.argmin(scores)
    min_rmsd_idx = np.argmin(ligand)

    df = pd.DataFrame({"min_score_name": [names[min_score_idx]], "min_score_soap": [scores[min_score_idx]], "min_score_ligand_rmsd": [ligand[min_score_idx]],
                       "min_score_interface_rmsd": [interface[min_score_idx]],
                       "min_rmsd_name": [names[min_rmsd_idx]], "min_rmsd_soap": [scores[min_rmsd_idx]], "min_rmsd_ligand_rmsd": [ligand[min_rmsd_idx]],
                       "min_rmsd_interface_rmsd": [interface[min_rmsd_idx]], "diff_rmsd":[ligand[min_score_idx] - ligand[min_rmsd_idx]]})

    df.to_csv("dock_summery_" + os.path.basename(directory) + ".csv", header=True, index=False)
    os.chdir("..")
    return df




if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="directory path containing the pdb directories")
    args = parser.parse_args()
    os.chdir(args.directory)
    i = 1
    first_pdb = True
    with open("dock_score_rmsd_summery.csv", 'w') as sum_file:
        for directory in os.listdir(args.directory):
            #  if the folder is pdb folder
            if os.path.isdir(directory) and re.fullmatch("[a-zA-Z0-9]{4}_[0-9]", directory):
                if directory in BAD_PDBS:  # TODO - delete when they are fixed
                    continue
                df = main(os.path.join(args.directory, directory))
                df["pdb_name"] = [directory]
                df.to_csv(sum_file, header=first_pdb, index=False)
                first_pdb = False
                print(i)
                i += 1










