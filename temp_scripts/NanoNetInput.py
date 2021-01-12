from Bio.PDB import *
import argparse
import os
import sys
import numpy as np
import pandas as pd
import pickle
from tqdm import tqdm
from NanoNetLabels import get_dist
from NanoNetUtils import *
import subprocess
import re

sys.path.insert(1, '/cs/usr/tomer.cohen13/lab/nanobodies/scripts')
from cdr_annotation import *


CDR_MAX_LENGTH = 32


def generate_input(pdb):
    """

    :param pdb:
    :return:
    """

    model_path = os.path.join(pdb, "ref.pdb")
    model = PDBParser().get_structure(id=pdb, file=model_path)[0]['H']

    seq, aa_residues = get_seq(model)
    cdr3_matrix = one_hot_coding(seq, 3)

    if "X" in seq:
        print("Warning, PDB: {}, has unknown aa".format(pdb))


    cdr1_matrix = one_hot_coding(seq, 1)

    # if OPTION == 1:
    #     third_matrix = option1(aa_residues, cdr1_start, cdr1_end, cdr3_start, cdr3_end)
    # elif OPTION == 2:  # OPTION == 2'
    #     pad = (CDR_MAX_LENGTH - (cdr1_end+1 - cdr1_start)) // 2
    #     third_matrix = option2(aa_residues, cdr1_start, cdr1_end, cdr3_start, cdr3_end, pad)
    # elif OPTION == 3:

    third_matrix = one_hot_coding(seq, 2)
    return np.dstack([cdr1_matrix, cdr3_matrix, third_matrix])


if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="dirctory path containing the pdb files")
    args = parser.parse_args()
    os.chdir(args.directory)
    failed_pdbs = pd.DataFrame(columns=["PDB", "FOLDER"])
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
                    input_matrix.append(generate_input(pdb))
            os.chdir("..")
    input_matrix = np.stack(input_matrix, axis=0)
    input_file_name = "nn_input"
    # if TEST:
    #     input_file_name += "_test"
    pickle.dump(input_matrix, open(input_file_name + ".pkl", "wb"))
    failed_pdbs.to_csv("nn_input_failed_pdbs.csv")


########################################################################################################################
#                                                                                                                      #
#                                       different options for 3rd matrix                                               #
#                                                                                                                      #
########################################################################################################################

# DIM = 2
# TEST = True
# OPTION = 3


def option2(residues, cdr1_start, cdr1_end, cdr3_start, cdr3_end, pad):
    """
    :return:
    """

    cdr1_residues = residues[cdr1_start:cdr1_end+1]
    dist = np.zeros((CDR_MAX_LENGTH, 21))
    start = residues[cdr3_start-1]["CA"].get_vector()
    end = residues[cdr3_end]["CA"].get_vector()
    for i in range(len(cdr1_residues)):
        for j in range(21):
            c1 = 'CB'
            if 'CB' not in cdr1_residues[i]:  # GLY
                c1 = 'CA'
            dist[i+pad][j] = (cdr1_residues[i][c1].get_vector() - (start + (end-start)**(j/(CDR_MAX_LENGTH-1)))).norm()
    print(np.array(dist))
    return np.array(dist)


def calc_dist(residues, cdr1_start, cdr1_end, cdr3_start, cdr3_end):
    """
    calculates the distances between the cdr2 C-alpha atoms and cdr3 start and end atoms
    :param pep: peptide object
    :param cdr2_start: int
    :param cdr2_end: int
    :param cdr3_start: int
    :param cdr3_end: int
    :return: dist array1, dist array2
    """

    dist1 = []
    dist2 = []
    cdr3_start_atom = residues[cdr3_start-1]["CA"]
    cdr3_end_atom = residues[cdr3_end ]["CA"]

    for residue in residues[cdr1_start:cdr1_end+1]:
        dist1.append(cdr3_start_atom - residue["CA"])
        dist2.append(cdr3_end_atom - residue["CA"])

    return dist1, dist2


def calc_angles(residues, cdr1_start, cdr1_end, cdr3_start, cdr3_end):
    """
    calculates the angles between the cdr2 C-alpha atoms and cdr3 start and end N, C-alpha atoms
    :param pep: peptide object
    :param cdr2_start: int
    :param cdr2_end: int
    :param cdr3_start: int
    :param cdr3_end: int
    :return: angle array1, angle array2
    """

    angle1 = []
    angle2 = []
    cdr3_start_resi = residues[cdr3_start -1]
    cdr3_end_resi = residues[cdr3_end ]

    for residue in residues[cdr1_start:cdr1_end+1]:
        angle1.append(calc_angle(cdr3_start_resi["N"].get_vector(), cdr3_start_resi["CA"].get_vector(), residue["CA"].get_vector()))
        angle2.append(calc_angle(cdr3_end_resi["N"].get_vector(), cdr3_end_resi["CA"].get_vector(), residue["CA"].get_vector()))

    return angle1, angle2


def option1(residues, cdr1_start, cdr1_end, cdr3_start, cdr3_end):
    """

    :param pep:
    :param cdr2_start:
    :param cdr2_end:
    :param cdr3_start:
    :param cdr3_end:
    :return:
    """

    cdr1_len = (cdr1_end+1 - cdr1_start)
    cdr1_pad = (CDR_MAX_LENGTH - cdr1_len) // 2

    dist_start, dist_end = calc_dist(residues, cdr1_start, cdr1_end, cdr3_start, cdr3_end)
    angle_start, angle_end = calc_angles(residues, cdr1_start, cdr1_end, cdr3_start, cdr3_end)

    dist_angle_matrix = np.zeros((CDR_MAX_LENGTH, 21))

    for i in range(0, 20, 4):
        dist_angle_matrix[cdr1_pad:cdr1_pad+cdr1_len, i] = dist_start
        dist_angle_matrix[cdr1_pad:cdr1_pad+cdr1_len, i+1] = dist_end
        dist_angle_matrix[cdr1_pad:cdr1_pad+cdr1_len, i+2] = angle_start
        dist_angle_matrix[cdr1_pad:cdr1_pad+cdr1_len, i+3] = angle_end

    return dist_angle_matrix
