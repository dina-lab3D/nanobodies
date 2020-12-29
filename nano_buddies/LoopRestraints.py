

from NanoNetInut import generate_input
from modeller import *
from modeller.automodel import *
from modeller import soap_loop
from modeller.scripts import complete_pdb
import argparse
import numpy as np
from NanoNetUtils import remove_pad

import Bio
from Bio import PDB
from Bio import SeqIO
from Bio.Blast import NCBIXML
import cdr_annotation

import subprocess
import os,sys,getopt

current_home = os.path.dirname(sys.argv[0])
print (current_home)

blast_home = "/cs/labs/dina/dina/software/ncbi-blast-2.8.1+/bin/"
rmsd_prog = "/cs/staff/dina/utils/rmsd"
get_pdb = "/cs/staff/dina/scripts/getPDB.pl"
get_pdb_chains = "/cs/staff/dina/scripts/getPDBChains.pl"
renumber = "/cs/staff/dina/utils/srcs/renumber/renumber"
rmsd_align = "/cs/staff/dina/scripts/alignRMSD.pl"
get_frag_chain = "/cs/staff/dina/utils/get_frag_chain.Linux"

loop_model_num = 200


def loop_model(seq, pdb_file, restraints_matrix):


    pdb_file = pdb_file
    [cdr1_s, cdr1_e] = cdr_annotation.find_cdr1(seq)
    [cdr3_s, cdr3_e] = cdr_annotation.find_cdr3(seq)
    restraints_matrix = remove_pad(seq, restraints_matrix, 3)

    if (cdr3_e - cdr1_s + 1) != restraints_matrix.shape[0] or restraints_matrix.shape[0] != restraints_matrix.shape[1]:
        print("cdr3 sequence error!!!")

    class MyLoop(loopmodel):
    # This routine picks the residues to be refined by loop modeling
        def select_loop_atoms(self):
            return selection(self.residue_range(cdr3_s+1, cdr3_e-2), self.residue_range(cdr1_s+1, cdr1_e-2)) # focus  TODO: cdr1 loop modeling

        def special_restraints(self, aln):

            rsr = self.restraints
            at = self.atoms

            #  distances
            distance_restraints = restraints_matrix[0]
            for i in range(distance_restraints.shape[0]):
                for j in range(distance_restraints.shape[1]):

                    atom_i = 'CB:'
                    atom_j = 'CB:'
                    if seq[i] == 'G':
                        atom_i = 'CA:'
                    if seq[i] == 'G':
                        atom_j = 'CA:'
                    rsr.add(forms.gaussian(group=physical.xy_distance,feature=features.distance(at[atom_i + str(i)], at[atom_j + str(j)]),mean=distance_restraints[i,j], stdev=1))

            #  omegas
            omega_restraints = np.arctan2(restraints_matrix[1,0], restraints_matrix[1,1])
            for i in range(omega_restraints.shape[0]):
                for j in range(omega_restraints.shape[1]):
                    if seq[i] == 'G' or seq[j] == 'G':
                        continue
                    rsr.add(forms.gaussian(group=physical.dihedral,feature=features.dihedral(at['CA:' + str(i)], at['CB:' + str(i)], at['CB:' + str(j)], at['CA:' + str(j)]),mean=omega_restraints[i,j], stdev=0.25))

            #  thethas
            thetha_restraints = np.arctan2(restraints_matrix[2,0], restraints_matrix[2,1])
            for i in range(thetha_restraints.shape[0]):
                for j in range(thetha_restraints.shape[1]):
                    if seq[i] == 'G' or seq[j] == 'G':
                        continue
                    rsr.add(forms.gaussian(group=physical.dihedral,feature=features.dihedral(at['N:' + str(i)], at['CA:' + str(i)], at['CB:' + str(i)], at['CB:' + str(j)]),mean=thetha_restraints[i,j], stdev=0.15))

            #  phis
            phis_restraints = np.arctan2(restraints_matrix[3,0], restraints_matrix[3,1])
            for i in range(phis_restraints.shape[0]):
                for j in range(phis_restraints.shape[1]):
                    if seq[i] == 'G' or seq[j] == 'G':
                        continue
                    rsr.add(forms.gaussian(group=physical.angle,feature=features.angle(at['CA:' + str(i)], at['CB:' + str(i)], at['CB:' + str(j)]),mean=phis_restraints[i,j], stdev=0.1))

    return MyLoop(env,inimodel=pdb_file,sequence='NANO')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", help="fasta file")
    parser.add_argument("pdb", help="pdb_model")
    parser.add_argument("NanoNet", help="trained NanoNet model")
    args = parser.parse_args()


    model_seq = get_sequence(args.fasta)

    [cdr1_start, cdr1_end] = cdr_annotation.find_cdr1(model_seq)  # TODO- fix the cdr1 bug
    [cdr3_start, cdr3_end] = cdr_annotation.find_cdr3(model_seq)

    restraints = args.NanoNet.predict(generate_input(args.pdb))

    m = loop_model(model_seq, args.pdb, restraints)
    m.loop.starting_model= 1           # index of the first loop model
    m.loop.ending_model = loop_model_num   # index of the last loop model
    m.loop.md_level = refine.slow      # loop refinement method; this yields
    # models quickly but of low quality;
    # use refine.slow for better models
    m.make()

    # renumber for calculating cdrs rmsd (same length of models...)
    if os.path.isfile("ref.pdb"):
        subprocess.run(renumber + " ref.pdb > ref_renumber.pdb", shell=True)

        subprocess.run(get_frag_chain + " ref_renumber.pdb H " + str(cdr1_start) + " " + str(cdr1_end) + " > ref_cdr1.pdb", shell=True)
        subprocess.run(get_frag_chain + " ref_renumber.pdb H " + str(cdr2_start) + " " + str(cdr2_end) + " > ref_cdr2.pdb", shell=True)
        subprocess.run(get_frag_chain + " ref_renumber.pdb H " + str(cdr3_start) + " " + str(cdr3_end) + " > ref_cdr3.pdb", shell=True)


    # score models
    f = open("network_scores.txt", "w")
    cdr_dist = open("cdrs_dist", "w")


    model_name = 'model_0.pdb'
    mdl = complete_pdb(env, model_name)
    s = selection(mdl)   # all atom selection

    dope_score = s.assess_dope(output='ENERGY_PROFILE NO_REPORT',
                               normalize_profile=True, smoothing_window=15)
    soap_score = s.assess(sp, output='ENERGY_PROFILE NO_REPORT',
                          normalize_profile=True, smoothing_window=15)
    rmsd = 0.0
    cdr1_rmsd = 0.0
    cdr2_rmsd = 0.0
    cdr3_rmsd = 0.0
    if os.path.isfile("ref.pdb"):

        cmd = rmsd_prog + " -t ref.pdb " + model_name + " | tail -n1 "
        rmsd_out = subprocess.check_output(cmd, shell=True)
        rmsd = float(rmsd_out.strip())

        # cdr 1,2,3 rmsd

        subprocess.run(rmsd_align + " ref.pdb" + " " + model_name, shell=True)  # align to get rmsd of cdr without cheating...

        subprocess.run(get_frag_chain + " " + model_name.replace(".pdb", "_tr.pdb") + " ' ' " + str(cdr1_start) + " " + str(cdr1_end) + " > temp_cdr1.pdb", shell=True)
        subprocess.run(get_frag_chain + " " + model_name.replace(".pdb", "_tr.pdb") + " ' ' " + str(cdr2_start) + " " + str(cdr2_end) + " > temp_cdr2.pdb", shell=True)
        subprocess.run(get_frag_chain + " " + model_name.replace(".pdb", "_tr.pdb") + " ' ' " + str(cdr3_start) + " " + str(cdr3_end) + " > temp_cdr3.pdb", shell=True)

        cdr1_rmsd = float(subprocess.check_output(rmsd_prog + " ref_cdr1.pdb temp_cdr1.pdb | tail -n1 ", shell=True).strip())
        cdr2_rmsd = float(subprocess.check_output(rmsd_prog + " ref_cdr2.pdb temp_cdr2.pdb | tail -n1 ", shell=True).strip())
        cdr3_rmsd = float(subprocess.check_output(rmsd_prog + " ref_cdr3.pdb temp_cdr3.pdb | tail -n1 ", shell=True).strip())

        calc_dist(model_name.replace(".pdb", "_tr.pdb"), cdr3_start, cdr3_end, cdr_dist)
        os.remove(model_name.replace(".pdb", "_tr.pdb"))


    print ("MODEL ", model_name, " dope-score: ", dope_score, " soap-score: ", soap_score, " rmsd: ", rmsd, " cdr1-rmsd: ", 0, " cdr2-rmsd: ", cdr2_rmsd, " cdr3-rmsd: ", cdr3_rmsd)
    ""
    f.write("MODEL "+ model_name + " dope-score: " + str(dope_score) + " soap-score: " + str(soap_score) + " rmsd: " + str(rmsd) +
            " cdr1-rmsd: " + str(cdr1_rmsd) + " cdr2-rmsd: " + str(cdr2_rmsd) + " cdr3-rmsd: " + str(cdr3_rmsd) + "\n")


    # score loops
    for i in range(1, (loop_model_num)+1):  # change to range(1, (loop_model_num*2)+1) for loop modeling 2
        # read model file
        code = "NANO.BL%04d0001.pdb" % i
        mdl = complete_pdb(env, code)
        s = selection(mdl)
        dope_score = s.assess_dope(output='ENERGY_PROFILE NO_REPORT',
                                   normalize_profile=True, smoothing_window=15)
        soap_score = s.assess(sp, output='ENERGY_PROFILE NO_REPORT',
                              normalize_profile=True, smoothing_window=15)
        rmsd = 0.0
        cdr1_rmsd = 0.0
        cdr2_rmsd = 0.0
        cdr3_rmsd = 0.0

        if os.path.isfile("ref.pdb"):
            cmd = rmsd_prog + " -t ref.pdb " + code + " | tail -n1 "
            rmsd_out = subprocess.check_output(cmd, shell=True)
            rmsd = float(rmsd_out.strip())

        # cdr 1,2,3 rmsd

        subprocess.run(rmsd_align + " ref.pdb" + " " + code, shell=True)  # align to get rmsd of cdr without cheating...

        subprocess.run(get_frag_chain + " " + code.replace(".pdb", "_tr.pdb") + " ' ' " + str(cdr1_start) + " " + str(cdr1_end) + " > temp_cdr1.pdb", shell=True)
        subprocess.run(get_frag_chain + " " + code.replace(".pdb", "_tr.pdb") + " ' ' " + str(cdr2_start) + " " + str(cdr2_end) + " > temp_cdr2.pdb", shell=True)
        subprocess.run(get_frag_chain + " " + code.replace(".pdb", "_tr.pdb") + " ' ' " + str(cdr3_start) + " " + str(cdr3_end) + " > temp_cdr3.pdb", shell=True)

        cdr1_rmsd = float(subprocess.check_output(rmsd_prog + " ref_cdr1.pdb temp_cdr1.pdb | tail -n1 ", shell=True).strip())
        cdr2_rmsd = float(subprocess.check_output(rmsd_prog + " ref_cdr2.pdb temp_cdr2.pdb | tail -n1 ", shell=True).strip())
        cdr3_rmsd = float(subprocess.check_output(rmsd_prog + " ref_cdr3.pdb temp_cdr3.pdb | tail -n1 ", shell=True).strip())

        print ("LOOP ", code, " dope-score: ", dope_score, " soap-score: ", soap_score, " rmsd: ", rmsd, " cdr1-rmsd: ", 0, " cdr2-rmsd: ", cdr2_rmsd, " cdr3-rmsd: ", cdr3_rmsd)
        f.write("LOOP "+ code + " dope-score: " + str(dope_score) + " soap-score: " + str(soap_score) + " rmsd: " + str(rmsd) +
                " cdr1-rmsd: " + str(cdr1_rmsd) + " cdr2-rmsd: " + str(cdr2_rmsd) + " cdr3-rmsd: " + str(cdr3_rmsd) + "\n")

        file_to_remove = "NANO.DL%04d0001" % i
        os.remove(file_to_remove)
        os.remove(code.replace(".pdb", "_tr.pdb"))

    # clean cdrs files
    if os.path.isfile("ref.pdb"):
        os.remove("temp_cdr1.pdb")
        os.remove("temp_cdr2.pdb")
        os.remove("temp_cdr3.pdb")
        os.remove("ref_renumber.pdb")

    # # clean up templates
    # for template in template_list:
    #     code = template[0]
    #     chain = template[1]
    #     pdb = code + '.pdb'
    #     pdb_chain = code + chain + '.pdb'
    #     if os.path.exists(pdb):
    #         os.remove(pdb)
    #         os.remove(pdb_chain)
    #
    # # clean up other files
    # os.remove("NANO.rsr")
    # os.remove("NANO.ini")
    # os.remove("NANO.lrsr")
    # os.remove("NANO.sch")
    # os.remove("default")
