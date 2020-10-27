
import argparse
import os
import subprocess
import re


# max memory for each batch
MEMORY = "8000m"

# max time for each batch
TIME = "8:0:0"

# paths to scripts and programs that are used in PyDock
BUILD_PARAM = "/cs/staff/dina/projects2/PatchDock/buildParams.pl  "  # TODO - change if want params with restrections (to /cs/usr/tomer.cohen13/nanobodies/COVID_19/PatchDock/buildParams.pl)
GET_CHAIN = "~dina/utils/getChain.Linux "
PATCH_DOCK = "/cs/staff/dina/projects2/PatchDock/patch_dock.Linux "
PATCH_DOCK_TRANS = "/cs/staff/dina/projects2/PatchDock/PatchDockOut2Trans.pl "
SETUP_ENV = "/cs/labs/dina/dina/libs/imp_build/setup_environment.sh "
SOAP_SCORE = "/cs/labs/dina/dina/libs/imp_build/bin/soap_score "
DOCK_DATA_MAKER = "/cs/labs/dina/tomer.cohen13/nanobodies/nano_buddies/DockDataMaker.py "
CLUSTER = "/cs/labs/dina/tomer.cohen13/InterfaceClustering/interface_cluster.linux "
BEST_TRANS = "/cs/labs/dina/tomer.cohen13/nanobodies/COVID_19/CovidBestTrans.py "

# the begining of the script for cluster
INTRO = "#!/bin/tcsh\n" \
        "#SBATCH --mem="+ MEMORY +"\n" \
        "#SBATCH -c1\n" \
        "#SBATCH --time=" + TIME + "\n"

SPIKE_PDB = "/cs/usr/tomer.cohen13/lab/nanobodies/COVID_19/S1.pdb "  # TODO - change if want different Antigen


def intro(script_file):
    """
    writes the intro of the script file
    :param script_file: open script file
    :return: None
    """
    script_file.write(INTRO)
    script_file.write("cd " + os.getcwd() + "\n")
    script_file.write("module load opencv\nsetenv CGAL_DIR /cs/labs/dina/dina/libs/CGAL\n")


def dock_pdb(folder):
    """
    creates a script and runs the Patch_Dock algorithm on all the loop/model pdb files in the directory
    :param folder: pdb directory (after running NanobodySelector.py on the folder)
    :return: None
    """
    os.chdir(folder)
    pdb_folder = os.getcwd()
    pdb_name = os.path.basename(pdb_folder)

    #  build and save the nanobody pdb
    nanobody = pdb_name + "_nanobody.pdb"
    subprocess.run(GET_CHAIN + " H " + pdb_name + ".pdb" + " > " + nanobody, shell=True)

    #  write docking script for all the loop.pdb/model.pdb
    with open("dock_script.sh", 'w') as script_file:
        intro(script_file)

        #  build parameters list
        script_file.write(BUILD_PARAM + nanobody + " " + SPIKE_PDB + " 4.0 AA\n")

        #  change name parameters list
        params_name = "params_" + pdb_name + ".txt"
        script_file.write("mv params.txt " + params_name + "\n")

        script_file.write("cp mycdrs3 cdrs3" + "\n")
        script_file.write("cp myframe frame" + "\n")

        # Docking
        docking_name = "docking_" + pdb_name + ".res"
        script_file.write(PATCH_DOCK + params_name + " " + docking_name + "\n")
        trans_name = "trans_" + pdb_name
        script_file.write(PATCH_DOCK_TRANS + docking_name + " > " + trans_name + "\n")

        #  soap scores
        script_file.write(SETUP_ENV + SOAP_SCORE + SPIKE_PDB + " " + nanobody + " " + trans_name + " -o soap_score_" + pdb_name + ".res\n")

        # make dock_data.csv file
        script_file.write("python3 " + DOCK_DATA_MAKER + " " + os.getcwd() + "\n")
        # make soap_score_cluster.res file (clustering)
        script_file.write(CLUSTER + " -f " + SPIKE_PDB + " " + nanobody + " dock_data.csv 4 soap_score_cluster.res\n")

        # make file with the best transformations according to soap_score
        script_file.write("python3 " + BEST_TRANS + pdb_folder + "\n")

#  send the script to the cluster
    subprocess.run("sbatch dock_script.sh", shell=True)
    os.chdir("..")


if __name__ == '__main__':
    """
    runs the dock_pdb() function on all the pdbs folders in the given directory. make sure to run from cluster
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="directory path containing the pdb directories")
    args = parser.parse_args()
    os.chdir(args.directory)

    for directory in os.listdir(args.directory):
        #  if the folder is pdb folder
        if os.path.isdir(directory) and re.fullmatch("[a-zA-Z0-9]{4}_[0-9]", directory):
            dock_pdb(os.path.join(args.directory, directory))

