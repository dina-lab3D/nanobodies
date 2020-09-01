
import argparse
import os
import subprocess
import re


# max memory for each batch
MEMORY = "10000m"

# max time for each batch
TIME = "8:0:0"

SCRIPT_PATH = "/cs/labs/dina/tomer.cohen13/nanobodies/nano_buddies/nanobodies_script.py"
BUILD_PARAM = "/cs/staff/dina/projects2/PatchDock/buildParams.pl "
CHAIN_SELECTOR = "~dina/scripts/chainSelector.pl "
GET_CHAIN = "~dina/utils/getChain.Linux "
PATCH_DOCK = "/cs/staff/dina/projects2/PatchDock/patch_dock.Linux "
PATCH_DOCK_TRANS = "/cs/staff/dina/projects2/PatchDock/PatchDockOut2Trans.pl "
SETUP_ENV = "/cs/labs/dina/dina/libs/imp_build/setup_environment.sh "
SOAP_SCORE = "/cs/labs/dina/dina/libs/imp_build/bin/soap_score "
RMSD_ALIGN = "/cs/staff/dina/scripts/alignRMSD.pl "
COMBINED_REF = "/cs/labs/dina/tomer.cohen13/nanobodies/nano_buddies/plots1000/combine_refs.py "

# the begining of the script for cluster
INTRO = "#!/bin/tcsh\n" \
        "#SBATCH --mem="+ MEMORY +"\n" \
        "#SBATCH -c1\n" \
        "#SBATCH --time=" + TIME + "\n"


def dock_pdb(directory):
    """
    creates a script and runs the Patch_Dock algorithm on all the loop/model pdb files in the directory
    :param directory: pdb directory (after running NanobodySelector.py on the folder)
    :return: None
    """
    os.chdir(directory)
    antigen_pdb = os.path.basename(directory) + ".pdb"

    #  gets a string containing the chains letters corresponding to the antigen chains (all but H)
    antigen_chains = subprocess.run(CHAIN_SELECTOR + antigen_pdb, shell=True,
                                    capture_output=True, universal_newlines=True).stdout.replace("H", "").replace(" ", "").replace("\n", "")
    #  build antigen pdb
    subprocess.run(GET_CHAIN + antigen_chains + " " + antigen_pdb + " > " + "antigen.pdb", shell=True)

    #  write docking script
    with open("dock_script.sh", 'w') as script_file:
        script_file.write(INTRO)
        script_file.write("cd " + os.getcwd() + "\n")
        script_file.write("module load opencv\nsetenv CGAL_DIR /cs/labs/dina/dina/libs/CGAL\n")
        for pdb_file in os.listdir(os.getcwd()):
            #  loop/ model nanobody pdb
            if (pdb_file.startswith("model") or pdb_file.startswith("loop")) and pdb_file.endswith(".pdb") and "tr" not in pdb_file:
                pdb_name = pdb_file.split(".")[0]

                #  align to ref.pdb to get correct rmsd
                script_file.write(RMSD_ALIGN + "ref.pdb " + pdb_file + "\n")
                tr_pdb_file = pdb_name + "_tr.pdb"

                #  parameters list
                script_file.write(BUILD_PARAM + tr_pdb_file + " antigen.pdb 4.0 AA\n")

                #  change name parameters list
                params_name = "params_" + pdb_name + ".txt"
                script_file.write("mv params.txt " + params_name + "\n")

                #  TODO
                script_file.write("cp mycdrs3 cdrs3" + "\n")
                #  TODO
                script_file.write("cp myframe frame" + "\n")

                # Docking
                if True:  #  TODO
                    docking_name = "docking_" + pdb_name + ".res"
                    script_file.write(PATCH_DOCK + params_name + " " + docking_name + "\n")
                    trans_name = "trans_" + pdb_name
                    script_file.write(PATCH_DOCK_TRANS + docking_name + " > " + trans_name + "\n")

                #  soap scores
                if True:  #  TODO
                    script_file.write(SETUP_ENV + SOAP_SCORE + "antigen.pdb " + tr_pdb_file + " " + trans_name + " -o soap_score_" + pdb_name + ".res\n")
                    script_file.write(SETUP_ENV + SOAP_SCORE + "antigen.pdb " + tr_pdb_file + " -o no_trans_soap_score_" + pdb_name + ".res\n")

    script_file.write(SETUP_ENV + SOAP_SCORE + "antigen.pdb ref.pdb -o soap_score_ref.res\n")
    #  send the script to the cluster
    script_file.write("python3 " + COMBINED_REF + " " + os.getcwd())
    subprocess.run("sbatch dock_script.sh", shell=True)
    os.chdir("..")


if __name__ == '__main__':
    """
    runs the dock_pdb() function on all the pdbs folders in the given directory.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="directory path containing the pdb directories")
    args = parser.parse_args()
    os.chdir(args.directory)

    for directory in os.listdir(args.directory):
        #  if the folder is pdb folder
        if os.path.isdir(directory) and re.fullmatch("[a-zA-Z0-9]{4}_[0-9]", directory):
            dock_pdb(os.path.join(args.directory, directory))

