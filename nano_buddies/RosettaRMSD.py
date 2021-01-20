
import argparse
import numpy as np
from NanoNetUtils import remove_pad, get_sequence, get_seq
import cdr_annotation

import subprocess
import os,sys,getopt
from Bio.PDB import *
import pandas as pd
import re

rmsd_prog = "/cs/staff/dina/utils/rmsd"
renumber = "/cs/staff/dina/utils/srcs/renumber/renumber"
rmsd_align = "/cs/staff/dina/scripts/alignRMSD.pl"
get_frag_chain = "/cs/staff/dina/utils/get_frag_chain.Linux"

N_LOOPS = 200
PROBLEM = ["1XGR_1"]


def calc_rmsds(pdb_dir):

    os.chdir(pdb_dir)
    model = PDBParser().get_structure("grafting/model-0.relaxed.pdb", "grafting/model-0.relaxed.pdb")[0]["H"]
    model_seq, _ = get_seq(model)

    ref_model = PDBParser().get_structure("ref.pdb", "ref.pdb")[0]["H"]
    ref_seq, _ = get_seq(ref_model)

    start, end = re.search(model_seq, ref_seq).span()
    subprocess.run(get_frag_chain + " ref.pdb H " + str(start) + " " + str(end) + " > ref_same_size.pdb", shell=True)

    [cdr1_start, cdr1_end] = cdr_annotation.find_cdr1(model_seq)
    [cdr2_start, cdr2_end] = cdr_annotation.find_cdr2(model_seq)
    [cdr3_start, cdr3_end] = cdr_annotation.find_cdr3(model_seq)

    if pdb_dir in PROBLEM:
        correct_seq = get_sequence(pdb_dir + ".fa")
        [cdr3_start, cdr3_end] = cdr_annotation.find_cdr3(correct_seq)
        print("PROBLEM")
        print(pdb_dir)
        print(cdr_annotation.find_cdr3(model_seq))
        print(cdr_annotation.find_cdr3(correct_seq))

    subprocess.run(get_frag_chain + " ref_same_size.pdb H " + str(cdr1_start) + " " + str(cdr1_end) + " > ref_cdr1.pdb", shell=True)
    subprocess.run(get_frag_chain + " ref_same_size.pdb H " + str(cdr2_start) + " " + str(cdr2_end) + " > ref_cdr2.pdb", shell=True)
    subprocess.run(get_frag_chain + " ref_same_size.pdb H " + str(cdr3_start) + " " + str(cdr3_end) + " > ref_cdr3.pdb", shell=True)

    all_rmsds = []
    cdr1_rmsds = []
    cdr2_rmsds = []
    cdr3_rmsds = []
    cdr3_lengths = []

    # score loops
    for i in range(1, N_LOOPS+1):  # change to range(1, (loop_model_num*2)+1) for loop modeling 2
        # read model file
        code = "model-0.relaxed_%04d.pdb" % i
        model_path = os.path.join("H3_modeling", code)
        cmd = rmsd_prog + " -t ref_same_size.pdb " + model_path + " | tail -n1 "
        rmsd_out = subprocess.check_output(cmd, shell=True)
        rmsd = float(rmsd_out.strip())

        # cdr 1,2,3 rmsd
        subprocess.run(rmsd_align + " ref_same_size.pdb " + model_path, shell=True)  # align to get rmsd of cdr without cheating...
        subprocess.run(get_frag_chain + " " + model_path.replace(".pdb", "_tr.pdb") + " H " + str(cdr1_start) + " " + str(cdr1_end) + " > temp_cdr1.pdb", shell=True)
        subprocess.run(get_frag_chain + " " + model_path.replace(".pdb", "_tr.pdb") + " H " + str(cdr2_start) + " " + str(cdr2_end) + " > temp_cdr2.pdb", shell=True)
        subprocess.run(get_frag_chain + " " + model_path.replace(".pdb", "_tr.pdb") + " H " + str(cdr3_start) + " " + str(cdr3_end) + " > temp_cdr3.pdb", shell=True)

        cdr1_rmsd = float(subprocess.check_output(rmsd_prog + " ref_cdr1.pdb temp_cdr1.pdb | tail -n1 ", shell=True).strip())
        cdr2_rmsd = float(subprocess.check_output(rmsd_prog + " ref_cdr2.pdb temp_cdr2.pdb | tail -n1 ", shell=True).strip())
        cdr3_rmsd = float(subprocess.check_output(rmsd_prog + " ref_cdr3.pdb temp_cdr3.pdb | tail -n1 ", shell=True).strip())

        all_rmsds.append(rmsd)
        cdr1_rmsds.append(cdr1_rmsd)
        cdr2_rmsds.append(cdr2_rmsd)
        cdr3_rmsds.append(cdr3_rmsd)
        cdr3_lengths.append(cdr3_end - cdr3_start)
        os.remove(model_path.replace(".pdb", "_tr.pdb"))
        os.remove("temp_cdr1.pdb")
        os.remove("temp_cdr2.pdb")
        os.remove("temp_cdr3.pdb")

    rosetta_scores = pd.read_csv("H3_modeling_scores.fasc", skiprows=1, delim_whitespace=True)
    rosetta_scores.drop_duplicates(subset="description", keep="last", inplace=True) # when running rosetta in parallel it can cause several models with same name

    rosetta_scores["rmsd"] = all_rmsds
    rosetta_scores["cdr1_rmsd"] = cdr1_rmsds
    rosetta_scores["cdr2_rmsd"] = cdr2_rmsds
    rosetta_scores["cdr3_rmsd"] = cdr3_rmsds
    rosetta_scores["cdr3_length"] = cdr3_lengths
    rosetta_scores.rename(columns={"description":"name"}, inplace=True)
    rosetta_scores["type"] = ["ROSETTA_LOOP"] * len(rosetta_scores["rmsd"])
    rosetta_scores.to_csv("H3_modeling_scores_rmsd.csv", index=False, columns=["type", "name", "total_score", "unconstr_score", "score", "rmsd", "cdr1_rmsd", "cdr2_rmsd", "cdr3_rmsd", "cdr3_length"])

    os.chdir("..")


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("pdbs", help="pdsb_dir")
    args = parser.parse_args()

    os.chdir(args.pdbs)
    for pdb in os.listdir(os.getcwd()):
        if os.path.exists(os.path.join(pdb, "H3_modeling_scores.fasc")):
            calc_rmsds(pdb)
            print("{} ended successfully".format(pdb))
        else:
            print("{} Failed, no score file".format(pdb))
