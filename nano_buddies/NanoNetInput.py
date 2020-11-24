from Bio.PDB import *
import argparse
import os
import sys
import numpy as np
import pickle
from tqdm import tqdm
sys.path.insert(1, '/cs/usr/tomer.cohen13/lab/nanobodies/scripts')
from cdr_annotation import *

CDR_MAX_LENGTH = 32
AA_DICT = {"A": 0, "C": 1, "D": 2, "E": 3, "F": 4, "G": 5, "H": 6, "I": 7, "K": 8, "L": 9, "M": 10, "N": 11, "P": 12,
           "Q": 13, "R": 14, "S": 15, "T": 16, "W": 17, "Y": 18, "V": 19, "-": 20}

DIM = 2


def calc_dist(pep, cdr2_start, cdr2_end, cdr3_start, cdr3_end):
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
    cdr3_start_atom = pep[cdr3_start -1]["CA"]
    cdr3_end_atom = pep[cdr3_end ]["CA"]

    for residue in pep[cdr2_start:cdr2_end+1]:
        dist1.append(cdr3_start_atom - residue["CA"])
        dist2.append(cdr3_end_atom - residue["CA"])

    return dist1, dist2


def calc_angles(pep, cdr2_start, cdr2_end, cdr3_start, cdr3_end):
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
    cdr3_start_resi = pep[cdr3_start -1]
    cdr3_end_resi = pep[cdr3_end ]

    for residue in pep[cdr2_start:cdr2_end+1]:
        angle1.append(calc_angle(cdr3_start_resi["N"].get_vector(), cdr3_start_resi["CA"].get_vector(), residue["CA"].get_vector()))
        angle2.append(calc_angle(cdr3_end_resi["N"].get_vector(), cdr3_end_resi["CA"].get_vector(), residue["CA"].get_vector()))

    return angle1, angle2


def get_dist_angle_matrix(pep, cdr2_start, cdr2_end, cdr3_start, cdr3_end):
    """

    :param pep:
    :param cdr2_start:
    :param cdr2_end:
    :param cdr3_start:
    :param cdr3_end:
    :return:
    """

    cdr2_len = (cdr2_end+1 - cdr2_start)
    cdr2_pad = (CDR_MAX_LENGTH - cdr2_len) // 2

    dist_start, dist_end = calc_dist(pep, cdr2_start, cdr2_end, cdr3_start, cdr3_end)
    angle_start, angle_end = calc_angles(pep, cdr2_start, cdr2_end, cdr3_start, cdr3_end)

    dist_angle_matrix = np.zeros((CDR_MAX_LENGTH, 21))

    dist_angle_matrix[cdr2_pad:cdr2_pad+cdr2_len, 9] = dist_start
    dist_angle_matrix[cdr2_pad:cdr2_pad+cdr2_len, 10] = dist_end
    dist_angle_matrix[cdr2_pad:cdr2_pad+cdr2_len, 11] = angle_start
    dist_angle_matrix[cdr2_pad:cdr2_pad+cdr2_len, 12] = angle_end

    return dist_angle_matrix


def generate_input(pdb):
    """

    :param pdb:
    :return:
    """
    if not os.path.exists(os.path.join(pdb, "model_0.pdb")):
        return None
    model = PDBParser().get_structure(pdb, os.path.join(pdb, "model_0.pdb"))[0]
    pep = PPBuilder().build_peptides(model)[0]
    seq = str(pep.get_sequence())

    [cdr2_start, cdr2_end], [cdr3_start, cdr3_end] = find_cdr2(seq), find_cdr3(seq)

    cdr2_len, cdr3_len = (cdr2_end+1 - cdr2_start), (cdr3_end+1 - cdr3_start)

    cdr2_pad, cdr3_pad = (CDR_MAX_LENGTH - cdr2_len) // 2, (CDR_MAX_LENGTH - cdr3_len) // 2

    seq_cdr2 = cdr2_pad * "-" + seq[cdr2_start:cdr2_end+1] + (CDR_MAX_LENGTH - cdr2_pad - cdr2_len) * "-"
    seq_cdr3 = cdr3_pad * "-" + seq[cdr3_start:cdr3_end+1] + (CDR_MAX_LENGTH - cdr3_pad - cdr3_len) * "-"

    cdr2_matrix, cdr3_matrix = np.zeros((CDR_MAX_LENGTH, 21)), np.zeros((CDR_MAX_LENGTH, 21))

    for i in range(CDR_MAX_LENGTH):
        cdr2_matrix[i][AA_DICT[seq_cdr2[i]]] = 1
        cdr3_matrix[i][AA_DICT[seq_cdr3[i]]] = 1
    if DIM == 1:
        return cdr3_matrix
    dist_angle_matrix = get_dist_angle_matrix(pep, cdr2_start, cdr2_end, cdr3_start, cdr3_end)
    return np.stack([cdr2_matrix, cdr3_matrix, dist_angle_matrix], axis=0)


if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="dirctory path containing the pdb files")
    args = parser.parse_args()
    os.chdir(args.directory)

    input_matrix = []
    for directory in os.listdir(os.getcwd()):
        if os.path.isdir(directory) and int(directory) <=5:  # directories 1,2,3,4...
            os.chdir(directory)
            for pdb in tqdm(os.listdir(os.getcwd())):
                if os.path.isdir(pdb):
                    pdb_input = generate_input(pdb)
                    if pdb_input is None:
                        print(directory + ": " + pdb + ", FAILED")
                        continue
                    input_matrix.append(pdb_input)
            os.chdir("..")
    input_matrix = np.stack(input_matrix, axis=0)
    pickle.dump(input_matrix, open("nn_input_" + str(DIM) + ".pkl", "wb"))
