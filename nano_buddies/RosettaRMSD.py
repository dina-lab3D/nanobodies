
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
MATCH_SUB = []
MATCH_FULL = []
MATCH_CUT = []
SEQS = []


def adjust_lengths(ref_seq, model_seq, pdb_dir):

    subprocess.run(renumber + " ref.pdb > ref_renumber.pdb", shell=True, stdout=subprocess.DEVNULL)
    match = re.search(model_seq, ref_seq)
    indexes = None

    if match:
        start, end = match.span()
        subprocess.run(get_frag_chain + " ref_renumber.pdb H " + str(start+1) + " " + str(end) + " > temp.pdb", shell=True, stdout=subprocess.DEVNULL)
        subprocess.run(renumber + " temp.pdb > ref_same_size.pdb", shell=True, stdout=subprocess.DEVNULL)
        os.remove("temp.pdb")

    else:
        match = re.search(ref_seq, model_seq)
        if match:
            subprocess.run("cp ref_renumber.pdb ref_same_size.pdb", shell=True, stdout=subprocess.DEVNULL)
            indexes = match.span()

            MATCH_FULL.append(pdb_dir)
        else:
            match = re.search(ref_seq[:-1], model_seq)
            if match:
                subprocess.run(get_frag_chain + " ref_renumber.pdb H " + str(1) + " " + str(len(ref_seq) - 1) + " > ref_same_size.pdb", shell=True, stdout=subprocess.DEVNULL)
                indexes = match.span()

                MATCH_SUB.append(pdb_dir)
            else:
                i = len(model_seq)
                while not match:
                    i -= 1
                    match = re.search(model_seq[:i], ref_seq)
                    if i == "50":
                        print(pdb_dir)
                        exit(1)

                start, end = match.span()
                subprocess.run(get_frag_chain + " ref_renumber.pdb H " + str(start+1) + " " + str(end) + " > temp.pdb", shell=True, stdout=subprocess.DEVNULL)
                subprocess.run(renumber + " temp.pdb > ref_same_size.pdb", shell=True, stdout=subprocess.DEVNULL)
                os.remove("temp.pdb")
                indexes = (0, i)
                MATCH_CUT.append(pdb_dir)
                SEQS.append([ref_seq[start:end], model_seq[0:i]])

    os.remove("ref_renumber.pdb")
    return indexes


def calc_rmsds(pdb_dir, nano_net):

    os.chdir(pdb_dir)

    model = PDBParser().get_structure("H3_NanoNet_modeling/model-0.relaxed_0001.pdb", "H3_NanoNet_modeling/model-0.relaxed_0001.pdb")[0]["H"]
    model_seq, _ = get_seq(model)

    ref_model = PDBParser().get_structure("ref.pdb", "ref.pdb")[0]["H"]
    ref_seq, _ = get_seq(ref_model)
    indexes = adjust_lengths(ref_seq, model_seq, pdb_dir)

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

    if indexes:
        cdr1_start -= indexes[0]
        cdr2_start -= indexes[0]
        cdr3_start -= indexes[0]
        cdr1_end -= indexes[0]
        cdr2_end -= indexes[0]
        cdr3_end -= indexes[0]

    subprocess.run(get_frag_chain + " ref_same_size.pdb H " + str(cdr1_start+1) + " " + str(cdr1_end) + " > ref_cdr1.pdb", shell=True, stdout=subprocess.DEVNULL)
    subprocess.run(get_frag_chain + " ref_same_size.pdb H " + str(cdr2_start+1) + " " + str(cdr2_end) + " > ref_cdr2.pdb", shell=True, stdout=subprocess.DEVNULL)
    subprocess.run(get_frag_chain + " ref_same_size.pdb H " + str(cdr3_start+1) + " " + str(cdr3_end) + " > ref_cdr3.pdb", shell=True, stdout=subprocess.DEVNULL)

    all_rmsds = []
    cdr1_rmsds = []
    cdr2_rmsds = []
    cdr3_rmsds = []
    cdr3_lengths = []

    if nano_net:
        folder = "H3_NanoNet_modeling"
        score_file = "H3_NanoNet_modeling_scores.fasc"
        rmsd_file = "H3_NanoNet_modeling_scores_rmsd.csv"

    else:
        folder = "H3_modeling"
        score_file = "H3_modeling_scores.fasc"
        rmsd_file = "H3_modeling_scores_rmsd.csv"

    # score loops
    for i in range(1, N_LOOPS+1):  # change to range(1, (loop_model_num*2)+1) for loop modeling 2
        # read model file
        code = "model-0.relaxed_%04d.pdb" % i
        model_path = os.path.join(folder, code)
        subprocess.run(renumber + " {} > model_renumber.pdb".format(model_path), shell=True, stdout=subprocess.DEVNULL)
        if indexes:
            subprocess.run(get_frag_chain + " model_renumber.pdb H " + str(indexes[0]+1) + " " + str(indexes[1]) + " > temp.pdb", shell=True, stdout=subprocess.DEVNULL)
            os.remove("model_renumber.pdb")
            subprocess.run(renumber + " temp.pdb > model_renumber.pdb", shell=True, stdout=subprocess.DEVNULL)
            os.remove("temp.pdb")

        cmd = rmsd_prog + " -t ref_same_size.pdb model_renumber.pdb | tail -n1 "
        rmsd = float(subprocess.run(cmd, shell=True, capture_output=True).stdout.strip())

        # cdr 1,2,3 rmsd
        subprocess.run(rmsd_align + " ref_same_size.pdb model_renumber.pdb", shell=True, stdout=subprocess.DEVNULL)  # align to get rmsd of cdr without cheating...
        subprocess.run(get_frag_chain + " model_renumber_tr.pdb H " + str(cdr1_start+1) + " " + str(cdr1_end) + " > temp_cdr1.pdb", shell=True, stdout=subprocess.DEVNULL)
        subprocess.run(get_frag_chain + " model_renumber_tr.pdb H " + str(cdr2_start+1) + " " + str(cdr2_end) + " > temp_cdr2.pdb", shell=True, stdout=subprocess.DEVNULL)
        subprocess.run(get_frag_chain + " model_renumber_tr.pdb H " + str(cdr3_start+1) + " " + str(cdr3_end) + " > temp_cdr3.pdb", shell=True, stdout=subprocess.DEVNULL)

        cdr1_rmsd = float(subprocess.run(rmsd_prog + " ref_cdr1.pdb temp_cdr1.pdb | tail -n1 ", shell=True, capture_output=True).stdout.strip())
        cdr2_rmsd = float(subprocess.run(rmsd_prog + " ref_cdr2.pdb temp_cdr2.pdb | tail -n1 ", shell=True, capture_output=True).stdout.strip())
        cdr3_rmsd = float(subprocess.run(rmsd_prog + " ref_cdr3.pdb temp_cdr3.pdb | tail -n1 ", shell=True, capture_output=True).stdout.strip())

        all_rmsds.append(rmsd)
        cdr1_rmsds.append(cdr1_rmsd)
        cdr2_rmsds.append(cdr2_rmsd)
        cdr3_rmsds.append(cdr3_rmsd)
        cdr3_lengths.append(cdr3_end - cdr3_start)
        os.remove("model_renumber.pdb")
        os.remove("model_renumber_tr.pdb")
        os.remove("temp_cdr1.pdb")
        os.remove("temp_cdr2.pdb")
        os.remove("temp_cdr3.pdb")

    rosetta_scores = pd.read_csv(score_file, skiprows=1, delim_whitespace=True)
    rosetta_scores.drop_duplicates(subset="description", keep="last", inplace=True) # when running rosetta in parallel it can cause several models with same name
    # a = np.array(rosetta_scores["description"])
    # for i in range(1, N_LOOPS+1):  # change to range(1, (loop_model_num*2)+1) for loop modeling 2
    #     code = "model-0.relaxed_%04d" % i
    #     if code not in a:
    #         print(code)
    rosetta_scores["rmsd"] = all_rmsds
    rosetta_scores["cdr1_rmsd"] = cdr1_rmsds
    rosetta_scores["cdr2_rmsd"] = cdr2_rmsds
    rosetta_scores["cdr3_rmsd"] = cdr3_rmsds
    rosetta_scores["length"] = cdr3_lengths
    rosetta_scores.rename(columns={"description":"name"}, inplace=True)
    rosetta_scores["type"] = ["LOOP"] * len(rosetta_scores["rmsd"])
    rosetta_scores.to_csv(rmsd_file, index=False, columns=["type", "name", "total_score", "unconstr_score", "score", "rmsd", "cdr1_rmsd", "cdr2_rmsd", "cdr3_rmsd", "length"])

    os.chdir("..")


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("pdbs", help="pdsb_dir")
    parser.add_argument("-n", "--nano_net", help="use results after using restraints with NanoNet", action="store_true")
    args = parser.parse_args()

    k = 1
    os.chdir(args.pdbs)
    for pdb in os.listdir(os.getcwd()):

        if pdb != "2FL5_1":
            continue
        # os.chdir(pdb)
        # subprocess.run("rm -f ref_cdr*", shell=True)
        # subprocess.run("rm -f H3_modeling_scores_rmsd.csv", shell=True)
        # subprocess.run("rm -f ref_same_size.pdb", shell=True)
        # subprocess.run("rm -f model_renumber.pdb", shell=True)
        # subprocess.run("rm -f temp.pdb", shell=True)
        # subprocess.run("rm -f temp_cdr1.pdb", shell=True)
        # subprocess.run("rm -f temp_cdr2.pdb", shell=True)
        # subprocess.run("rm -f temp_cdr3.pdb", shell=True)
        #
        # os.chdir("..")
        if os.path.exists(os.path.join(pdb, "H3_modeling_scores.fasc")):
            print(pdb)
            calc_rmsds(pdb, args.nano_net)
            print("{} ended successfully".format(pdb))
        else:
            print("{} Failed, no score file".format(pdb))
        print("NUM: " + str(k))
        k += 1

    print(MATCH_FULL)
    print(MATCH_SUB)
    print(MATCH_CUT)
    print(SEQS)
