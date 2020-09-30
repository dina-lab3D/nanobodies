
import argparse
import os
import subprocess
import re
import pandas as pd
import numpy as np

# max memory for each batch
MEMORY = "10000m"

# max time for each batch
TIME = "8:0:0"

BUILD_PARAM = "/cs/staff/dina/projects2/PatchDock/buildParams.pl "
BUILD_PARAM_XL = "/cs/staff/dina/projects2/PatchDock/buildParamsXlinksAA.pl "
CHAIN_SELECTOR = "~dina/scripts/chainSelector.pl "
GET_CHAIN = "~dina/utils/getChain.Linux "
PATCH_DOCK = "/cs/staff/dina/projects2/PatchDock/patch_dock.Linux "
PATCH_DOCK_TRANS = "/cs/staff/dina/projects2/PatchDock/PatchDockOut2Trans.pl "
SETUP_ENV = "/cs/labs/dina/dina/libs/imp_build/setup_environment.sh "
SOAP_SCORE = "/cs/labs/dina/dina/libs/imp_build/bin/soap_score "
RMSD_ALIGN = "/cs/staff/dina/scripts/alignRMSD.pl "
DOCK_DATA_MAKER = "/cs/labs/dina/tomer.cohen13/nanobodies/nano_buddies/DockDataMaker.py "
CLUSTER = "/cs/labs/dina/tomer.cohen13/InterfaceClustering/interface_cluster.linux "
SET_CHAIN_NAME = "~dina/scripts/namechain.pl "

# the begining of the script for cluster
INTRO = "#!/bin/tcsh\n" \
        "#SBATCH --mem="+ MEMORY +"\n" \
        "#SBATCH -c1\n" \
        "#SBATCH --time=" + TIME + "\n"

NOISE = 0.25


def cross_links_selector(pdb_name, xl_num, noise):

    xl_file = os.path.join(os.path.dirname(os.getcwd()), "Jwalk_results", pdb_name + "_crosslink_list.txt")
    df = pd.read_csv(xl_file, sep="\s+", header=None, skiprows=[0], names=["Index", "Model", "Atom1", "Atom2", "SASD", "Euclidean Distance"])
    df["SASD"] = np.asfarray(df["SASD"], float)
    df["Euclidean Distance"] = np.asfarray(df["Euclidean Distance"], float)
    df["Index"] = np.asfarray(df["Index"], int)

    good_xl = df[df["SASD"] <= 25]
    bad_xl = df[df["SASD"] > 25]

    good_xl_indexs = np.array(good_xl["Index"]) - 1
    bad_xl_indexes = np.array(bad_xl["Index"]) - 1
    np.random.shuffle(good_xl_indexs), np.random.shuffle(bad_xl_indexes)

    noise_num = int(xl_num * noise)
    xl_df = pd.concat([good_xl.loc[good_xl_indexs[:xl_num-noise_num], ["Atom1", "Atom2"]], bad_xl.loc[bad_xl_indexes[:noise_num], ["Atom1", "Atom2"]]], ignore_index=True)

    xl_final = pd.DataFrame({"atom1_n": xl_df.Atom1.replace({"[A-Z]{3}-(\d+)-[A-Z]-[A-Z]+": r"\1"}, regex=True),
                             "atom1": xl_df.Atom1.replace({"[A-Z]{3}-\d+-([A-Z])-[A-Z]+": r"\1"}, regex=True),
                             "atom2_n": xl_df.Atom2.replace({"[A-Z]{3}-(\d+)-[A-Z]-[A-Z]+": r"\1"}, regex=True),
                             "atom2": xl_df.Atom2.replace({"[A-Z]{3}-\d+-([A-Z])-[A-Z]+": r"\1"}, regex=True),
                             "low": [0] * len(xl_df.Atom2), "high": [25] * len(xl_df.Atom2)})
    xl_final.to_csv("dist_constraints", sep=" ", header=False, index=False)


def dock_pdb(directory, xl):
    """
    creates a script and runs the Patch_Dock algorithm on all the loop/model pdb files in the directory
    :param directory: pdb directory (after running NanobodySelector.py on the folder)
    :return: None
    """
    os.chdir(directory)
    pdb_folder = os.getcwd()
    antigen_pdb = os.path.basename(directory) + ".pdb"

    #  gets a string containing the chains letters corresponding to the antigen chains (all but H)
    antigen_chains = subprocess.run(CHAIN_SELECTOR + antigen_pdb, shell=True,
                                    capture_output=True, universal_newlines=True).stdout.replace("H", "").replace(" ", "").replace("\n", "")
    #  build antigen pdb
    subprocess.run(GET_CHAIN + antigen_chains + " " + antigen_pdb + " > " + "antigen.pdb", shell=True)

    docking_folder = "pydock_results"
    if xl:
        docking_folder += "_xl_" + str(xl)
    if not os.path.exists(docking_folder):
        os.mkdir(docking_folder)
    os.chdir(docking_folder)

    #  write docking script
    with open("dock_script.sh", 'w') as script_file:
        script_file.write(INTRO)
        script_file.write("cd " + os.getcwd() + "\n")
        script_file.write("module load opencv\nsetenv CGAL_DIR /cs/labs/dina/dina/libs/CGAL\n")
        script_file.write("rm cat_soap_scores.res\n")  # remove cat_soap_scores.res if already exists (so we can append later)
        script_file.write("echo -n > cat_soap_scores.res\n")
        # for pdb_file in os.listdir(pdb_folder):
        #     #  loop/ model nanobody pdb
        #     if (pdb_file.startswith("model") or pdb_file.startswith("loop")) and pdb_file.endswith(".pdb") and "tr" not in pdb_file:
        #         pdb_name = pdb_file.split(".")[0]
        #
        #         #  align to ref.pdb to get correct rmsd
        #         script_file.write(RMSD_ALIGN + os.path.join(pdb_folder, "ref.pdb") + " " + os.path.join(pdb_folder, pdb_file) + "\n")
        #         tr_pdb_file = os.path.join(pdb_folder, pdb_name + "_tr.pdb")
        #
        #         #  parameters list
        #         if xl:
        #             cross_links_selector(os.path.basename(directory), xl, noise=NOISE)
        #             script_file.write(SET_CHAIN_NAME + tr_pdb_file + " H\n")
        #             script_file.write(BUILD_PARAM_XL + tr_pdb_file + " " + os.path.join(pdb_folder, "antigen.pdb") + " " + str(0.5) + " \n")
        #
        #         else:
        #             script_file.write(BUILD_PARAM + tr_pdb_file + " " + os.path.join(pdb_folder, "antigen.pdb") +" 4.0 AA\n")
        #         #  change name parameters list
        #         params_name = "params_" + pdb_name + ".txt"
        #         script_file.write("mv params.txt " + params_name + "\n")
        #
        #         #  TODO
        #         script_file.write("cp mycdrs3 cdrs3" + "\n")
        #         #  TODO
        #         script_file.write("cp myframe frame" + "\n")
        #
        #         # Docking
        #         if True:  #  TODO
        #             docking_name = "docking_" + pdb_name + ".res"
        #             script_file.write(PATCH_DOCK + params_name + " " + docking_name + "\n")
        #             trans_name = "trans_" + pdb_name
        #             script_file.write(PATCH_DOCK_TRANS + docking_name + " > " + trans_name + "\n")
        #
        #         #  soap scores
        #         if True:  #  TODO
        #             script_file.write(SETUP_ENV + SOAP_SCORE + os.path.join(pdb_folder, "antigen.pdb") + " " + tr_pdb_file + " " + trans_name + " -o soap_score_" + pdb_name + ".res\n")
        #             script_file.write(SETUP_ENV + SOAP_SCORE + os.path.join(pdb_folder, "antigen.pdb") + " " + tr_pdb_file + " -o no_trans_soap_score_" + pdb_name + ".res\n")
        #             script_file.write("grep \"|\" soap_score_" + pdb_name + ".res | grep -v SOAP | cut -d '|' -f1,2,7 >> cat_soap_scores.res\n")

        ###############################################

        script_file.write(BUILD_PARAM + os.path.join(pdb_folder, "ref.pdb") + " " + os.path.join(pdb_folder, "antigen.pdb") +" 4.0 AA\n")
        params_name = "params_ref.txt"
        script_file.write("mv params.txt " + params_name + "\n")
        script_file.write("cp mycdrs3 cdrs3" + "\n")
        script_file.write("cp myframe frame" + "\n")
        # Docking
        if True:  #  TODO
            docking_name = "docking_ref.res"
            script_file.write(PATCH_DOCK + params_name + " " + docking_name + "\n")
            trans_name = "trans_ref"
            script_file.write(PATCH_DOCK_TRANS + docking_name + " > " + trans_name + "\n")
        #  soap scores
        if True:  #  TODO
            script_file.write(SETUP_ENV + SOAP_SCORE + os.path.join(pdb_folder, "antigen.pdb") + " " + os.path.join(pdb_folder, "ref.pdb") + " " + trans_name + " -o soap_score_ref.res\n")

    ######################################
        # script_file.write(SETUP_ENV + SOAP_SCORE + os.path.join(pdb_folder, "antigen.pdb") + " " + os.path.join(pdb_folder, "ref.pdb") + " -o soap_score_ref_temp.res\n")
        # flag = " "
        # if xl:
        #     flag = " -xl"
        # script_file.write("python3 " + DOCK_DATA_MAKER + " " + os.getcwd() + flag + "\n")
        # script_file.write(CLUSTER + " -f " + os.path.join(pdb_folder, "antigen.pdb") + " " + os.path.join(pdb_folder, "ref.pdb") + " dock_data.csv 4 soap_score_cluster.res\n")

    #  send the script to the cluster
    subprocess.run("sbatch dock_script.sh", shell=True)
    os.chdir("..")
    os.chdir("..")


if __name__ == '__main__':
    """
    runs the dock_pdb() function on all the pdbs folders in the given directory.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="directory path containing the pdb directories")
    parser.add_argument("-xl", "--crosslinks", help="use cross links for docking", type=int)
    args = parser.parse_args()
    os.chdir(args.directory)

    for directory in os.listdir(args.directory):
        #  if the folder is pdb folder
        if os.path.isdir(directory) and re.fullmatch("[a-zA-Z0-9]{4}_[0-9]", directory):
            dock_pdb(os.path.join(args.directory, directory), args.crosslinks)

