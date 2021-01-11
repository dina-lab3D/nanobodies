from tensorflow.keras.models import load_model
from NanoNetUtils import generate_input
import argparse
import numpy as np
from NanoNetUtils import remove_pad, get_sequence
import cdr_annotation
import os


DIST_STD = 0.7
OMEGA_STD = 2
THETA_STD = 1.7
PHI_STD = 0


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
            const_file.write("AtomPair {} {} {} {} HARMONIC {} {}\n".format(atom_i_type, atom_i, atom_j_type, atom_j, constraints[i,j], DIST_STD))


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
            const_file.write("Dihedral CA {} CB {} CB {} CA {} CIRCULARHARMONIC {} {}\n".format(atom_i, atom_i, atom_j, atom_j, constraints[i, j], OMEGA_STD))


def write_const_theta(const_file, constraints, cdr3_s, seq):

    length = len(constraints)
    for i in range(length):
        for j in range(length):
            if i == j:
                continue
            atom_i = i + cdr3_s
            atom_j = j + cdr3_s
            if seq[atom_i] == 'G' or seq[atom_j] == 'G':  # GLY
                continue
            atom_i += 1  # pdb numbering starts from 1
            atom_j += 1  # pdb numbering starts from 1
            const_file.write(
                "Dihedral N {} CA {} CB {} CB {} CIRCULARHARMONIC {} {}\n".format(atom_i, atom_i, atom_i, atom_j, constraints[i, j], THETA_STD))


def write_const_phi(const_file, constraints, cdr3_s, seq):
    length = len(constraints)
    for i in range(length):
        for j in range(length):
            if i == j:
                continue
            atom_i = i + cdr3_s
            atom_j = j + cdr3_s
            if seq[atom_i] == 'G' or seq[atom_j] == 'G':  # GLY
                continue
            atom_i += 1  # pdb numbering starts from 1
            atom_j += 1  # pdb numbering starts from 1
            const_file.write("Dihedral CA {} CB {} CB {} CIRCULARHARMONIC {} {}\n".format(atom_i, atom_i, atom_j, constraints[i, j], PHI_STD))


def sanity_check(distance_restraints, cdr3_s, cdr3_e):

    print(distance_restraints.shape)
    if (cdr3_e - cdr3_s + 1) != distance_restraints.shape[0] or distance_restraints.shape[0] != distance_restraints.shape[1]:
        print(cdr3_e - cdr3_s + 1)
        print(distance_restraints.shape[0])
        print(distance_restraints.shape[1])
        print("cdr3 sequence error!!!")
        exit(2)

    if np.sum(distance_restraints[0]) == 0 or np.sum(distance_restraints[-1]) == 0 \
            or np.sum(distance_restraints[:,0]) ==0 or np.sum(distance_restraints[:,-1]) == 0:
        print("remove pad error!!!")
        exit(3)


def write_const_file(fasta_file, restraints_matrix):

    seq = get_sequence(fasta_file)

    distance_restraints = remove_pad(restraints_matrix[0][0,:,:,0], seq)
    omega_restraints = np.arctan2(remove_pad(restraints_matrix[1][0, :, :, 0], seq),remove_pad(restraints_matrix[1][0, :, :, 1], seq))
    thetha_restraints = np.arctan2(remove_pad(restraints_matrix[2][0, :, :, 0], seq),remove_pad(restraints_matrix[2][0, :, :, 1], seq))
    phis_restraints = np.arctan2(remove_pad(restraints_matrix[3][0, :, :, 0], seq),remove_pad(restraints_matrix[3][0, :, :, 1], seq))

    cdr_s, cdr_e = cdr_annotation.find_cdr3(seq)
    sanity_check(distance_restraints, cdr_s, cdr_e)

    with open(pdb_dir + "constraints", 'w') as const_file:
        write_const_dist(const_file, distance_restraints, cdr_s, seq)
        write_const_omega(const_file, omega_restraints, cdr_s, seq)
        write_const_theta(const_file, thetha_restraints, cdr_s, seq)
        write_const_phi(const_file, phis_restraints, cdr_s, seq)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("Rosetta_pdbs", help="Rosetta PDBs directory")
    parser.add_argument("NanoNet", help="trained NanoNet model")
    args = parser.parse_args()

    nano_net_model = load_model(args.NanoNet)
    os.chdir(args.Rosetta_pdbs)

    for pdb_dir in os.listdir(os.getcwd()):
        os.chdir(pdb_dir)
        pdb_fasta = pdb_dir + ".fa"
        restraints = nano_net_model.predict(np.array([generate_input(pdb_fasta)]))
        write_const_file(pdb_fasta, restraints)
        os.chdir("..")


