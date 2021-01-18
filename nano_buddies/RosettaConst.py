from tensorflow.keras.models import load_model
from NanoNetUtils import generate_input
import argparse
import numpy as np
from NanoNetUtils import remove_pad, get_sequence, get_seq
import cdr_annotation
import os
from Bio.PDB import *

DIST_STD = 0.63
OMEGA_STD = 0.424
THETA_STD = 0.3
PHI_STD = 0.22

PROBLEM = ["1XGR_1"]


def write_const_dist(const_file, constraints, cdr3_s, seq):

    length = len(constraints)
    for i in range(length):
        for j in range(i+1, length):  # symmetry
            atom_i = i + cdr3_s
            atom_j = j + cdr3_s
            atom_i_type = 'CA' if seq[atom_i] == 'G' else 'CB'  # GLY
            atom_j_type = 'CA' if seq[atom_j] == 'G' else 'CB'  # GLY
            atom_i += 1  # pdb numbering starts from 1
            atom_j += 1  # pdb numbering starts from 1
            const_file.write("AtomPair {} {} {} {} HARMONIC {:.5f} {}\n".format(atom_i_type, atom_i, atom_j_type, atom_j, constraints[i,j], DIST_STD))


def write_const_omega(const_file, constraints, cdr3_s, seq):

    length = len(constraints)
    for i in range(length):
        for j in range(i + 1, length):  # symmetry
            atom_i = i + cdr3_s
            atom_j = j + cdr3_s
            if seq[atom_i] == 'G' or seq[atom_j] == 'G':  # GLY
                continue
            atom_i += 1  # pdb numbering starts from 1
            atom_j += 1  # pdb numbering starts from 1
            const_file.write("Dihedral CA {} CB {} CB {} CA {} CIRCULARHARMONIC {:.5f} {}\n".format(atom_i, atom_i, atom_j, atom_j, constraints[i, j], OMEGA_STD))


def write_const_theta(const_file, constraints, cdr3_s, seq):

    length = len(constraints)
    for i in range(length):
        for j in range(length):
            if i == j:  # same atom...
                continue
            atom_i = i + cdr3_s
            atom_j = j + cdr3_s
            if seq[atom_i] == 'G' or seq[atom_j] == 'G':  # GLY
                continue
            atom_i += 1  # pdb numbering starts from 1
            atom_j += 1  # pdb numbering starts from 1
            const_file.write("Dihedral N {} CA {} CB {} CB {} CIRCULARHARMONIC {:.5f} {}\n".format(atom_i, atom_i, atom_i, atom_j, constraints[i, j], THETA_STD))


def write_const_phi(const_file, constraints, cdr3_s, seq):
    length = len(constraints)
    for i in range(length):
        for j in range(length):
            if i == j:  # same atom...
                continue
            atom_i = i + cdr3_s
            atom_j = j + cdr3_s
            if seq[atom_i] == 'G' or seq[atom_j] == 'G':  # GLY
                continue
            atom_i += 1  # pdb numbering starts from 1
            atom_j += 1  # pdb numbering starts from 1
            const_file.write("Angle CA {} CB {} CB {} CIRCULARHARMONIC {:.5f} {}\n".format(atom_i, atom_i, atom_j, constraints[i, j], PHI_STD))


def write_const_file(sequence, restraints_matrix):

    distance_restraints = remove_pad(restraints_matrix[0][0,:,:,0], sequence) * 10  # we divided by factor 10 in NanoNet
    omega_restraints = np.arctan2(remove_pad(restraints_matrix[1][0, :, :, 1], sequence), remove_pad(restraints_matrix[1][0, :, :, 0], sequence))  # angle = arctan(sin, cos)
    thetha_restraints = np.arctan2(remove_pad(restraints_matrix[2][0, :, :, 1], sequence), remove_pad(restraints_matrix[2][0, :, :, 0], sequence))
    phis_restraints = np.arctan2(remove_pad(restraints_matrix[3][0, :, :, 1], sequence), remove_pad(restraints_matrix[3][0, :, :, 0], sequence))

    cdr_s, cdr_e = cdr_annotation.find_cdr3(sequence)
    # print(cdr_e - cdr_s+1)
    with open(pdb_dir + "_constraints", 'w') as const_file:
        write_const_dist(const_file, distance_restraints, cdr_s, sequence)
        write_const_omega(const_file, omega_restraints, cdr_s, sequence)
        write_const_theta(const_file, thetha_restraints, cdr_s, sequence)
        write_const_phi(const_file, phis_restraints, cdr_s, sequence)


def sanity_check(all_restraints, pdb, sequence):

    pdb_model = PDBParser().get_structure(pdb, "grafting/model-0.relaxed.pdb")[0]["H"]
    seq_fa = get_sequence(pdb + ".fa")
    seq_aa, _ = get_seq(pdb_model)
    cdr3_s, cdr3_e = cdr_annotation.find_cdr3(sequence)
    distance_restraints = remove_pad(all_restraints[0][0, :, :, 0], sequence) * 10
    error = False

    if (cdr3_e - cdr3_s + 1) != distance_restraints.shape[0] or distance_restraints.shape[0] != distance_restraints.shape[1]:
        print(cdr3_e - cdr3_s + 1)
        print(distance_restraints.shape[0])
        print(distance_restraints.shape[1])
        print("cdr3 sequence error!!!")
        error = True

    if np.sum(distance_restraints[0]) == 0 or np.sum(distance_restraints[-1]) == 0 \
            or np.sum(distance_restraints[:,0]) ==0 or np.sum(distance_restraints[:,-1]) == 0:
        print("remove pad error!!!")
        error = True

    a, b = cdr_annotation.find_cdr3(seq_aa)
    c, d = cdr_annotation.find_cdr3(seq_fa)

    if seq_aa[a:b+1] != seq_fa[c:d+1]:
        print(seq_aa[a:b+1])
        print(seq_fa[c:d+1])
        print(sequence[cdr3_s:cdr3_e+1])
        print("find_cdr3 error!!!")
        error = True
    if not error:
        print("{} finished successfully!".format(pdb))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("Rosetta_pdbs", help="Rosetta PDBs directory")
    parser.add_argument("NanoNet", help="trained NanoNet model")
    args = parser.parse_args()

    nano_net_model = load_model(args.NanoNet)
    os.chdir(args.Rosetta_pdbs)

    for pdb_dir in os.listdir(os.getcwd()):
        if os.path.isdir(pdb_dir) and pdb_dir != "1YC7_1":
            os.chdir(pdb_dir)
            model = PDBParser().get_structure(pdb_dir, "grafting/model-0.relaxed.pdb")[0]["H"]
            if pdb_dir in PROBLEM:
                seq = get_sequence(pdb_dir + ".fa")
                print(pdb_dir)
            else:
                seq, aa_chain = get_seq(model)

            restraints = nano_net_model.predict(np.array([generate_input(seq, fasta=False)]))
            write_const_file(seq, restraints)
            sanity_check(restraints, pdb_dir, seq)
            os.chdir("..")



