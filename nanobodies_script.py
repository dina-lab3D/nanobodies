import os
import subprocess
import sys


def main(abs_path_folder, file_name):

    folder_name = file_name.split(".")[0]
    # create nanobody folder
    if not os.path.isdir(abs_path_folder):
        os.mkdir(abs_path_folder)
    os.chdir(abs_path_folder)

    # get chain H (antibody)
    subprocess.run("~dina/utils/getChain.Linux H " + abs_path_folder + ".pdb" + " > " + abs_path_folder + "/ref.pdb", shell=True)

    # move nanobody to its folder
    # subprocess.run("mv " + os.path.abspath(file_name) + " " + folder_name + "/", shell=True)

    # fasta script
    subprocess.run("~dina/utils/pdb2fasta " + abs_path_folder + "/ref.pdb" + " > " + abs_path_folder + "/" + folder_name + ".fa",
                   shell=True)

    # run the script that creates the loops models
    os.chdir(abs_path_folder)

    subprocess.run("~dina/modeller9.18/bin/modpy.sh python /cs/labs/dina/tomer.cohen13/nanobodies/scripts/modelNanobody.py -t -l 1000" + abs_path_folder +
                   "/" + folder_name + ".fa", shell=True)




if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
