from Bio.PDB import *
import argparse
import os
import sys
import numpy as np
import pandas as pd
import pickle
from tqdm import tqdm
from NanoNetUtils import *

sys.path.insert(1, '/cs/usr/tomer.cohen13/lab/nanobodies/scripts')
from cdr_annotation import *


CDR_MAX_LENGTH = 32
TEST = True


def get_dist(pep_residues, start, end, pad=0):
    """

    :param pep_residues:
    :param start:
    :param end:
    :param pad:
    :return:
    """
    residues = pep_residues[start:end+1]
    dist = np.zeros((CDR_MAX_LENGTH, CDR_MAX_LENGTH))
    for i in range(len(residues)):
        for j in range(len(residues)):
            if i == j:
                continue
            c1 = 'CB'
            c2 = 'CB'
            if 'CB' not in residues[i]:  # GLY
                c1 = 'CA'
            if 'CB' not in residues[j]:  # GLY
                c2 = 'CA'
            dist[i+pad][j+pad] = (residues[i][c1] - residues[j][c2])
    return np.array([dist, dist])


def get_theta(pep_residues, start, end, pad=0):
    """

    :param pep_residues:
    :param start:
    :param end:
    :param pad:
    :return:
    """
    residues = pep_residues[start:end+1]
    cos_theta = np.zeros((CDR_MAX_LENGTH, CDR_MAX_LENGTH))
    sin_theta = np.zeros((CDR_MAX_LENGTH, CDR_MAX_LENGTH))

    for i in range(len(residues)):
        for j in range(len(residues)):
            if i == j:
                continue
            if not residues[i].has_id('CB') or not residues[j].has_id('CB'):  # GLY OR UNK OR MISSING CB
                continue

            angle = calc_dihedral(residues[i]["N"].get_vector(), residues[i]["CA"].get_vector(),
                                  residues[i]["CB"].get_vector(), residues[j]["CB"].get_vector())
            cos_theta[i+pad][j+pad] = np.cos(angle)
            sin_theta[i+pad][j+pad] = np.sin(angle)

    return np.array([cos_theta, sin_theta])


def get_phi(pep_residues, start, end, pad=0):
    """

    :param pep_residues:
    :param start:
    :param end:
    :param pad:
    :return:
    """
    residues = pep_residues[start:end+1]
    cos_phi = np.zeros((CDR_MAX_LENGTH, CDR_MAX_LENGTH))
    sin_phi = np.zeros((CDR_MAX_LENGTH, CDR_MAX_LENGTH))

    for i in range(len(residues)):
        for j in range(len(residues)):
            if i == j:
                continue
            if not residues[i].has_id('CB') or not residues[j].has_id('CB'):  # GLY OR UNK OR MISSING CB
                continue
            angle = calc_angle(residues[i]["CA"].get_vector(), residues[i]["CB"].get_vector(),
                               residues[j]["CB"].get_vector())

            cos_phi[i+pad][j+pad] = np.cos(angle)
            sin_phi[i+pad][j+pad] = np.sin(angle)
    return np.array([cos_phi, sin_phi])


def get_omega(pep_residues, start, end, pad=0):
    """

    :param pep_residues:
    :param start:
    :param end:
    :param pad:
    :return:
    """
    residues = pep_residues[start:end+1]
    cos_omega = np.zeros((CDR_MAX_LENGTH, CDR_MAX_LENGTH))
    sin_omega = np.zeros((CDR_MAX_LENGTH, CDR_MAX_LENGTH))

    for i in range(len(residues)):
        for j in range(len(residues)):
            if i == j:
                continue
            if not residues[i].has_id('CB') or not residues[j].has_id('CB'):  # GLY OR UNK OR MISSING CB
                continue
            angle = calc_dihedral(residues[i]["CA"].get_vector(), residues[i]["CB"].get_vector(),
                                  residues[j]["CB"].get_vector(), residues[j]["CA"].get_vector())
            cos_omega[i+pad][j+pad] = np.cos(angle)
            sin_omega[i+pad][j+pad] = np.sin(angle)

    return np.array([cos_omega, sin_omega])


def generate_label(pdb):
    """

    :param pdb:
    :return:
    """

    if not os.path.exists(os.path.join(pdb, "ref.pdb")):
        return None
    model = PDBParser().get_structure(pdb, os.path.join(pdb, "ref.pdb"))[0]["H"]
    # pep = PPBuilder().build_peptides(model, aa_only=False)[0]

    for i in model.get_residues():
        if not i.has_id("N"):
            print("no side chains")
            return None
        break

    seq, aa_residues = get_seq(model)


    [cdr3_start, cdr3_end] = find_cdr3(seq)

    # for padding the result matrix with zeros
    pad = (CDR_MAX_LENGTH - (cdr3_end+1 - cdr3_start)) // 2

    # get angles and distance
    theta = get_theta(aa_residues, cdr3_start, cdr3_end, pad)
    dist = get_dist(aa_residues, cdr3_start, cdr3_end, pad)
    phi= get_phi(aa_residues, cdr3_start, cdr3_end, pad)
    omega = get_omega(aa_residues, cdr3_start, cdr3_end, pad)

    if "X" in seq:
        print("Warning, PDB: {}, has unknown aa".format(pdb))

    return np.array([dist, omega, theta, phi])


if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="dirctory path containing the pdb files")
    args = parser.parse_args()
    os.chdir(args.directory)
    failed_pdbs = pd.DataFrame(columns=["PDB", "FOLDER"])
    feature_matrix = []
    for directory in os.listdir(os.getcwd()):
        if os.path.isdir(directory) and directory != "failed":  # directories 1,2,3,4...
            os.chdir(directory)
            print(directory)
            for pdb in tqdm(os.listdir(os.getcwd())):
                if os.path.isdir(pdb):
                    if not valid_pdb(pdb, TEST):
                        print(directory + ": " + pdb + ", FAILED")
                        failed_pdbs = failed_pdbs.append(pd.DataFrame({"PDB":[pdb], "FOLDER":[directory]}))
                        continue
                    feature_matrix.append(generate_label(pdb))
            os.chdir("..")
    feature_matrix = np.stack(feature_matrix, axis=0)
    labels_file_name = "nn_labels_"
    if TEST:
        labels_file_name += "_test"
    pickle.dump(feature_matrix, open(labels_file_name + ".pkl", "wb"))
    failed_pdbs.to_csv("nn_labels_failed_pdbs.csv")


# questions
# 1. padding? 2. angles (gly), 3.dist ca not cb

# empty line, zero all or one ?
