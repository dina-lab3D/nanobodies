from Bio.PDB import *
import os
import sys
import numpy as np
from Bio import SeqIO

sys.path.insert(1, '/cs/usr/tomer.cohen13/lab/nanobodies/scripts')
from cdr_annotation import *
# from modelNanobody import get_sequence

CDR_MAX_LENGTH = 32
AA_DICT = {"A": 0, "C": 1, "D": 2, "E": 3, "F": 4, "G": 5, "H": 6, "I": 7, "K": 8, "L": 9, "M": 10, "N": 11, "P": 12,
           "Q": 13, "R": 14, "S": 15, "T": 16, "W": 17, "Y": 18, "V": 19, "-": 20, "X": 20}


def get_sequence(fasta_filename):
    for seq_record in SeqIO.parse(fasta_filename, "fasta"):
        sequence = str(seq_record.seq)
        return sequence


def one_hot_coding(seq, cdr):
    """
    :param seq:
    :param cdr:
    :return:
    """
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


def valid_pdb(pdb, test=False):

    # if not os.path.exists(os.path.join(pdb, "ref.pdb")):
    #     return False
    model = PDBParser().get_structure(id=pdb, file=pdb)[0]["H"]
    for i in model.get_residues():
        if not i.has_id("N"):
            print(pdb + ": no side chains")
            return False
        break

    # if not test:
    #     if not os.path.exists(os.path.join(pdb, "model_0.pdb")):
    #         print(pdb + ": failed modeling")
    #         return False

    return True


def get_seq(chain):

    aa_residues = []
    seq = ""

    for residue in chain.get_residues():
        aa = residue.get_resname()
        if not is_aa(aa) or not residue.has_id('CA'):
            continue
        elif aa == "UNK":
            seq += "X"
            aa_residues.append(residue)
        else:
            seq += Polypeptide.three_to_one(residue.get_resname())
            aa_residues.append(residue)

    return seq, aa_residues


def remove_pad(one_hot_matrix, seq):

    [cdr_start, cdr_end] = find_cdr3(seq)

    cdr_len = (cdr_end + 1 - cdr_start)
    pad_left = (CDR_MAX_LENGTH - cdr_len) // 2
    pad_right = pad_left + cdr_len
    return one_hot_matrix[pad_left:pad_right, pad_left:pad_right]


def generate_input(pdb_fasta):
    """

    :param pdb:
    :return:
    """

    seq = get_sequence(pdb_fasta)
    cdr3_matrix = one_hot_coding(seq, 3)

    if "X" in seq:
        print("Warning, PDB: {}, has unknown aa".format(pdb_fasta))

    cdr1_matrix = one_hot_coding(seq, 1)
    third_matrix = one_hot_coding(seq, 2)

    return np.dstack([cdr1_matrix, cdr3_matrix, third_matrix])


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
    return np.dstack([sin_theta, cos_theta])


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
    return np.dstack([sin_phi, cos_phi])


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
    return np.dstack([sin_omega, cos_omega])


def generate_label(fasta, pdb):
    """

    :param fasta:
    :param pdb:
    :return:
    """

    model = PDBParser().get_structure(pdb, pdb)[0]["H"]
    seq, aa_residues = get_seq(model)

    seq_test = get_sequence(fasta)
    if seq_test != seq:
        print(seq)
        print(seq_test)
        exit(1)

    [cdr3_start, cdr3_end] = find_cdr3(seq)

    # for padding the result matrix with zeros
    pad = (CDR_MAX_LENGTH - (cdr3_end+1 - cdr3_start)) // 2

    # get angles and distance
    theta = get_theta(aa_residues, cdr3_start, cdr3_end, pad)
    dist = get_dist(aa_residues, cdr3_start, cdr3_end, pad)
    phi= get_phi(aa_residues, cdr3_start, cdr3_end, pad)
    omega = get_omega(aa_residues, cdr3_start, cdr3_end, pad)

    if "X" in seq:
        print("Warning, PDB: {}, has unknown aa".format(fasta))
    labels_matrix = np.array([dist, omega, theta, phi])
    # if BINS:
    #     labels_matrix = matrix_to_bins(labels_matrix)
    return labels_matrix

