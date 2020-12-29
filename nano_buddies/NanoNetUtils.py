from Bio.PDB import *
import os
import sys
import numpy as np

sys.path.insert(1, '/cs/usr/tomer.cohen13/lab/nanobodies/scripts')
from cdr_annotation import *
from modelNanobody import get_sequence

CDR_MAX_LENGTH = 32
AA_DICT = {"A": 0, "C": 1, "D": 2, "E": 3, "F": 4, "G": 5, "H": 6, "I": 7, "K": 8, "L": 9, "M": 10, "N": 11, "P": 12,
           "Q": 13, "R": 14, "S": 15, "T": 16, "W": 17, "Y": 18, "V": 19, "-": 20, "X": 20}


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

    if not os.path.exists(os.path.join(pdb, "ref.pdb")):
        return False
    model = PDBParser().get_structure(id=pdb, file=os.path.join(pdb, "ref.pdb"))[0]["H"]
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


def remove_pad(seq, one_hot_matrix, cdr):

    find = [find_cdr1, find_cdr2, find_cdr3]
    [cdr_start, cdr_end] = find[cdr - 1](seq)

    cdr_len = (cdr_end + 1 - cdr_start)

    pad_left = (CDR_MAX_LENGTH - cdr_len) // 2
    pad_right = pad_left - cdr_len

    return one_hot_matrix[pad_left:pad_right, pad_left:pad_right]