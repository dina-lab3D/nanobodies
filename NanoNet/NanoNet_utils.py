from Bio.PDB import *
import os
import sys
import numpy as np
from Bio import SeqIO

from cdr_annotation import *
CDR_MAX_LENGTH = 32
AA_DICT = {"A": 0, "C": 1, "D": 2, "E": 3, "F": 4, "G": 5, "H": 6, "I": 7, "K": 8, "L": 9, "M": 10, "N": 11, "P": 12,
           "Q": 13, "R": 14, "S": 15, "T": 16, "W": 17, "Y": 18, "V": 19, "-": 20, "X": 20}


def get_sequence(fasta_filename):
    for seq_record in SeqIO.parse(fasta_filename, "fasta"):
        sequence = str(seq_record.seq)
        return sequence


def one_hot_coding(seq, cdr):

    find = [find_cdr1, find_cdr2, find_cdr3]
    [cdr_start, cdr_end] = find[cdr - 1](seq)

    cdr_len = (cdr_end + 1 - cdr_start)

    cdr_pad = (CDR_MAX_LENGTH - cdr_len) // 2

    seq_cdr = cdr_pad * "-" + seq[cdr_start:cdr_end + 1] + (
                CDR_MAX_LENGTH - cdr_pad - cdr_len) * "-"

    cdr_matrix = np.zeros((CDR_MAX_LENGTH, 21))

    for i in range(CDR_MAX_LENGTH):
        cdr_matrix[i][AA_DICT[seq_cdr[i]]] = 1

    return cdr_matrix


def remove_pad(one_hot_matrix, seq):

    [cdr_start, cdr_end] = find_cdr3(seq)

    cdr_len = (cdr_end + 1 - cdr_start)
    pad_left = (CDR_MAX_LENGTH - cdr_len) // 2
    pad_right = pad_left + cdr_len
    return one_hot_matrix[pad_left:pad_right, pad_left:pad_right]


def generate_input(seq):

    cdr3_matrix = one_hot_coding(seq, 3)
    if "X" in seq:
        print("Warning, PDB has unknown aa")

    cdr1_matrix = one_hot_coding(seq, 1)
    third_matrix = one_hot_coding(seq, 2)

    return np.dstack([cdr1_matrix, cdr3_matrix, third_matrix])


def normalize_dist(dist):
    dist = np.clip(dist,0 , 20)
    return dist / 10


def get_dist(pep_residues, start, end, pad=0):

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

    dist = normalize_dist(dist)
    return np.dstack([dist, dist])


def get_theta(pep_residues, start, end, pad=0):

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

    return np.dstack([cos_theta, sin_theta])


def get_phi(pep_residues, start, end, pad=0):

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

    return np.dstack([cos_phi, sin_phi])


def get_omega(pep_residues, start, end, pad=0):

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

    return np.dstack([cos_omega, sin_omega])


def generate_label(pdb, cdr=3):

    model = PDBParser().get_structure(pdb, pdb)[0]["H"]
    seq, aa_residues = get_seq(model)

    find = [find_cdr1, find_cdr2, find_cdr3]
    [cdr_start, cdr_end] = find[cdr-1](seq)

    # for padding the result matrix with zeros
    pad = (CDR_MAX_LENGTH - (cdr_end+1 - cdr_start)) // 2

    # get angles and distance
    theta = get_theta(aa_residues, cdr_start, cdr_end, pad)
    dist = get_dist(aa_residues, cdr_start, cdr_end, pad)
    phi= get_phi(aa_residues, cdr_start, cdr_end, pad)
    omega = get_omega(aa_residues, cdr_start, cdr_end, pad)

    if "X" in seq:
        print("Warning, PDB: {}, has unknown aa".format(pdb))
    labels_matrix = np.array([dist, omega, theta, phi])

    return labels_matrix


