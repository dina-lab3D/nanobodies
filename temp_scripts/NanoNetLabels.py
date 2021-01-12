from Bio.PDB import *
import argparse
import os
import sys
import numpy as np
import pandas as pd
import pickle
from tqdm import tqdm
import re
from NanoNetUtils import *

sys.path.insert(1, '/cs/usr/tomer.cohen13/lab/nanobodies/scripts')
from cdr_annotation import *

CDR_MAX_LENGTH = 32


def normalize_dist(dist):
    dist = np.clip(dist,0,20)
    return dist / 10


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
    # if BINS:
    #     return dist
    dist = normalize_dist(dist)
    return np.dstack([dist, dist])


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

    angles = np.zeros((CDR_MAX_LENGTH, CDR_MAX_LENGTH))

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
            angles[i + pad][j + pad] = np.degrees(angle) % 360
    # if BINS:
    #     return angles
    return np.dstack([cos_theta, sin_theta])


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
    angles = np.zeros((CDR_MAX_LENGTH, CDR_MAX_LENGTH))

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
            angles[i+pad][j+pad] = np.degrees(angle) % 360
    # if BINS:
    #     return angles
    return np.dstack([cos_phi, sin_phi])


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
    angles = np.zeros((CDR_MAX_LENGTH, CDR_MAX_LENGTH))

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
            angles[i+pad][j+pad] = np.degrees(angle) % 360
    # if BINS:
    #     return angles
    return np.dstack([cos_omega, sin_omega])


def generate_label(pdb):
    """

    :param pdb:
    :return:
    """

    model = PDBParser().get_structure(pdb, os.path.join(pdb, "ref.pdb"))[0]["H"]
    # pep = PPBuilder().build_peptides(model, aa_only=False)[0]

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
    labels_matrix = np.array([dist, omega, theta, phi])
    # if BINS:
    #     labels_matrix = matrix_to_bins(labels_matrix)
    return labels_matrix


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
                    pdb_names.append(os.path.join(directory,pdb))
            os.chdir("..")
    feature_matrix = np.stack(feature_matrix, axis=0)
    labels_file_name = "nn_labels"
    pdb_names_file = "pdb_names"
    # if TEST:
    #     labels_file_name += "_test"
    #     pdb_names_file += "_test"
    # if BINS:
    #     labels_file_name += "_bins"
    #     pdb_names_file += "pdb_names"

    pickle.dump(feature_matrix, open(labels_file_name + ".pkl", "wb"))
    pickle.dump(np.array(pdb_names), open(pdb_names_file + ".pkl", "wb"))
    failed_pdbs.to_csv("nn_labels_failed_pdbs.csv")

########################################################################################################################
#                                                                                                                      #
#                                          for bins representation and                                                 #
#                                                                                                                      #
########################################################################################################################


# TEST = True
# BINS = False


def distance_bins():
    cur_bin = 0
    bin_width = 0.5
    bins = [(cur_bin+ bin_width*i, cur_bin+bin_width + bin_width*i) for i in range(40)]
    bins.append((20, np.inf))
    return bins


def omega_theta_bins():
    cur_bin = 0
    bin_width = 15
    bins = [(cur_bin + bin_width * i, cur_bin + bin_width + bin_width * i) for i in range(24)]
    return bins


def phi_bins():
    cur_bin = 0
    bin_width = 15
    bins = [(cur_bin + bin_width * i, cur_bin + bin_width + bin_width * i) for i in range(12)]
    return bins


def matrix_to_bins(labels_matrix):

    bins_matrix = np.zeros_like(labels_matrix, dtype=np.int)
    all_bins = [distance_bins(), omega_theta_bins(), omega_theta_bins(), phi_bins()]

    for feature in range(len(all_bins)):
        bins = all_bins[feature]
        for j in range(len(bins)):
            (low, high) = bins[j]
            mask = np.logical_and(labels_matrix[feature] >= low, labels_matrix[feature] < high)
            bins_matrix[feature][mask] = j
    return bins_matrix

