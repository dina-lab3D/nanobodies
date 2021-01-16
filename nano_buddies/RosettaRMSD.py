
import argparse
import numpy as np
from NanoNetUtils import remove_pad, get_sequence
import cdr_annotation

import subprocess
import os,sys,getopt


rmsd_prog = "/cs/staff/dina/utils/rmsd"
renumber = "/cs/staff/dina/utils/srcs/renumber/renumber"
rmsd_align = "/cs/staff/dina/scripts/alignRMSD.pl"
get_frag_chain = "/cs/staff/dina/utils/get_frag_chain.Linux"

N_LOOPS = 200


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="pdb_dir")
    args = parser.parse_args()

    os.chdir(args.pdb)

    model_seq = get_sequence(os.path.basename(args.pdb) + ".fa")

    [cdr1_start, cdr1_end] = cdr_annotation.find_cdr1(model_seq)  # TODO- fix the cdr1 bug
    [cdr3_start, cdr3_end] = cdr_annotation.find_cdr3(model_seq)
    [cdr2_start, cdr2_end] = cdr_annotation.find_cdr2(model_seq)

    # renumber for calculating cdrs rmsd (same length of models...)
    subprocess.run(renumber + " ref.pdb > ref_renumber.pdb", shell=True)

    all_rmsds = []
    cdr1_rmsds = []
    cdr2_rmsds = []
    cdr3_rmsds = []

    # score loops
    for i in range(1, N_LOOPS+1):  # change to range(1, (loop_model_num*2)+1) for loop modeling 2
        # read model file
        code = "model-0.relaxed_%04d.pdb" % i
        path = os.path.join("H3_modeling", code)
        cmd = rmsd_prog + " -t ref.pdb " + code + " | tail -n1 "
        rmsd_out = subprocess.check_output(cmd, shell=True)
        rmsd = float(rmsd_out.strip())

        # cdr 1,2,3 rmsd
        subprocess.run(rmsd_align + " ref.pdb " + code, shell=True)  # align to get rmsd of cdr without cheating...
        subprocess.run(get_frag_chain + " " + code.replace(".pdb", "_tr.pdb") + " ' ' " + str(cdr1_start) + " " + str(cdr1_end) + " > temp_cdr1.pdb", shell=True)
        subprocess.run(get_frag_chain + " " + code.replace(".pdb", "_tr.pdb") + " ' ' " + str(cdr2_start) + " " + str(cdr2_end) + " > temp_cdr2.pdb", shell=True)
        subprocess.run(get_frag_chain + " " + code.replace(".pdb", "_tr.pdb") + " ' ' " + str(cdr3_start) + " " + str(cdr3_end) + " > temp_cdr3.pdb", shell=True)

        cdr1_rmsd = float(subprocess.check_output(rmsd_prog + " ref_cdr1.pdb temp_cdr1.pdb | tail -n1 ", shell=True).strip())
        cdr2_rmsd = float(subprocess.check_output(rmsd_prog + " ref_cdr2.pdb temp_cdr2.pdb | tail -n1 ", shell=True).strip())
        cdr3_rmsd = float(subprocess.check_output(rmsd_prog + " ref_cdr3.pdb temp_cdr3.pdb | tail -n1 ", shell=True).strip())

        all_rmsds.append(rmsd)
        cdr1_rmsds.append(cdr1_rmsd)
        cdr2_rmsds.append(cdr2_rmsd)
        cdr3_rmsds.append(cdr3_rmsd)

        os.remove(code.replace(".pdb", "_tr.pdb"))
        os.remove("temp_cdr1.pdb")
        os.remove("temp_cdr2.pdb")
        os.remove("temp_cdr3.pdb")

    rosetta_scores = pd.read_csv("H3_modeling_scores.fasc", skiprows=1, delim_whitespace=True)
    rosetta_scores.drop([0,1]) # first two loops missing cause running rosetta parallel

    rosetta_scores["rmsd"] = all_rmsds
    rosetta_scores["cdr1_rmsd"] = cdr1_rmsds
    rosetta_scores["cdr2_rmsd"] = cdr2_rmsds
    rosetta_scores["cdr3_rmsd"] = cdr3_rmsds

    rosetta_scores.to_csv("H3_modeling_scores_rmsd.csv", index=False)
    print("ended successfully")
