import os
import subprocess
import sys

# -t for testing mode, -l <int> for number of loops to model
FLAGS = "-t -l 100"


def main(abs_path_folder, file_name):
    """
    runs the modelNanobody.py program on the file_name (pdb file) after 
    creating fasta file for chain H. all the files are saved in the folder
    path abs_path_folder
    :param abs_path_folder: the folder path of the file_name 
    :param file_name: the name of the pdb file (nanobody)
    :return: None
    """
    folder_name = file_name.split(".")[0]

    # get chain H (antibody)
    subprocess.run("~dina/utils/getChain.Linux H " + abs_path_folder + ".pdb" + " > " + abs_path_folder + "/ref.pdb", shell=True)

    # fasta script
    subprocess.run("~dina/utils/pdb2fasta " + abs_path_folder + "/ref.pdb" + " > " + abs_path_folder + "/" + folder_name + ".fa",
                   shell=True)

    # run the script that creates the loops models
    subprocess.run("~dina/modeller9.18/bin/modpy.sh python3 /cs/labs/dina/tomer.cohen13/nanobodies/scripts/modelNanobody.py " + FLAGS + " " + abs_path_folder +
                   "/" + folder_name + ".fa", shell=True)

    # move the pdb file
    subprocess.run("mv " + abs_path_folder + ".pdb " + abs_path_folder, shell=True)


if __name__ == '__main__':
    # sys.argv[1] = folder path
    # sys.argv[2] = pdb file name
    main(sys.argv[1], sys.argv[2])
