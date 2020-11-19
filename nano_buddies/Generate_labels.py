from Bio.PDB import *
import argparse
import os
import sys
import numpy as np
import pickle

sys.path.insert(1, '/cs/usr/tomer.cohen13/lab/nanobodies/scripts')
from cdr_annotation import *

CDR_MAX_LENGTH = 30


def get_dist(pep, start, end, pad=0):
    """

    :param pep:
    :param start:
    :param end:
    :param pad:
    :return:
    """
    residues = pep[start:end+1]
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
    return dist


def get_theta(pep, start, end, pad=0):
    """

    :param pep:
    :param start:
    :param end:
    :param pad:
    :return:
    """
    residues = pep[start:end+1]
    theta = np.zeros((CDR_MAX_LENGTH, CDR_MAX_LENGTH))

    for i in range(len(residues)):
        for j in range(len(residues)):
            if i == j:
                continue
            if residues[i].get_resname() == 'GLY':   # TODO: make zero?
                continue
            atom = "CB"
            if residues[j].get_resname() == 'GLY':
                atom = "CA"
            theta[i+pad][j+pad] = calc_dihedral(residues[i]["N"].get_vector(), residues[i]["CA"].get_vector(),
                                                residues[i]["CB"].get_vector(), residues[j][atom].get_vector())
    return theta


def get_phi(pep, start, end, pad=0):
    """

    :param pep:
    :param start:
    :param end:
    :param pad:
    :return:
    """
    residues = pep[start:end+1]
    phi = np.zeros((CDR_MAX_LENGTH, CDR_MAX_LENGTH))

    for i in range(len(residues)):
        for j in range(len(residues)):
            if i == j:
                continue
            if residues[i].get_resname() == 'GLY':   # TODO: make zero?
                continue
            atom = "CB"
            if residues[j].get_resname() == 'GLY':
                atom = "CA"
            phi[i+pad][j+pad] = calc_angle(residues[i]["CA"].get_vector(), residues[i]["CB"].get_vector(),
                                           residues[j][atom].get_vector())
    return phi


def get_omega(pep, start, end, pad=0):
    """

    :param pep:
    :param start:
    :param end:
    :param pad:
    :return:
    """
    residues = pep[start:end+1]
    omega = np.zeros((CDR_MAX_LENGTH, CDR_MAX_LENGTH))

    for i in range(len(residues)):
        for j in range(len(residues)):
            if i == j:
                continue
            if residues[i].get_resname() == 'GLY' or residues[j].get_resname() == 'GLY':    # TODO: make zero?
                continue
            omega[i+pad][j+pad] = calc_dihedral(residues[i]["CA"].get_vector(), residues[i]["CB"].get_vector(),
                                                residues[j]["CB"].get_vector(), residues[j]["CA"].get_vector())
    return omega


def generate_label(pdb):
    """

    :param pdb:
    :return:
    """
    model = PDBParser().get_structure(pdb, pdb)[0]
    pep = PPBuilder().build_peptides(model)[0]
    seq = str(pep.get_sequence())

    [cdr3_start, cdr3_end] = find_cdr3(seq)

    # for padding the result matrix with zeros
    pad = (CDR_MAX_LENGTH - (cdr3_end+1 - cdr3_start)) // 2

    # get angles and distance
    theta = get_theta(pep, cdr3_start, cdr3_end, pad)
    dist = get_dist(pep, cdr3_start, cdr3_end, pad)
    phi= get_phi(pep, cdr3_start, cdr3_end, pad)
    omega = get_omega(pep, cdr3_start, cdr3_end, pad)

    return np.stack([dist, omega, theta, phi], axis=0)


if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="dirctory path containing the pdb files")
    args = parser.parse_args()
    os.chdir(args.directory)

    feature_matrix = []
    for pdb in os.listdir(os.getcwd()):
        feature_matrix.append(generate_label(pdb))

    feature_matrix = np.stack(feature_matrix, axis=0)
    pickle.dump(feature_matrix, open("features", "wb"))

# questions
# 1. padding? 2. angles (gly), 3.dist ca not cb

# empty line, zero all or one ?
