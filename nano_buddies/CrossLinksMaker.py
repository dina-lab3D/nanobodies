import argparse
import os
import subprocess
import re

# max memory for each batch
MEMORY = "10000m"

# max time for each batch
TIME = "8:0:0"

# the begining of the script for cluster
INTRO = "#!/bin/tcsh\n" \
        "#SBATCH --mem="+ MEMORY +"\n" \
        "#SBATCH -c1\n" \
        "#SBATCH --time=" + TIME + "\n"

# path to the xlm_tools (checks if the cross links generated in INTERFACE_XL are realistic
XLM_TOOLS = "/cs/labs/dina/tomer.cohen13/XLM-Tools-master/xlmtools.v1.0.py "

# path to interface_cross_links - makes cross links list
INTERFACE_XL = "/cs/staff/dina/libs/imp_build/bin/interface_cross_links "


def make_cross_links(pdf_folder):
    """

    :param pdf_folder:
    :return:
    """
    os.chdir(pdf_folder)
    pdb_name = os.path.basename(pdf_folder) + ".pdb"
    with open("xl_script.sh", 'w') as script_file:
        # intro
        script_file.write(INTRO)
        script_file.write("cd " + os.getcwd() + "\n")
        script_file.write("module load opencv\nsetenv CGAL_DIR /cs/labs/dina/dina/libs/CGAL\n")

        # make cross links
        script_file.write(INTERFACE_XL + " antigen.pdb ref.pdb 25 -n -e\n")

        # validates cross links
        script_file.write("python3 " + XLM_TOOLS + "-xl_list " + os.path.join(pdf_folder, "cxms_all.dat") + " -pdb " + pdb_name + " -jwalk -depth\n")

    # send to cluster
    subprocess.run("sbatch xl_script.sh", shell=True)
    os.chdir("..")


if __name__ == '__main__':

    """
    runs the make_cross_links() function on all the pdbs folders in the given directory. make sure to run from cluster
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="directory path containing the pdb directories")
    args = parser.parse_args()
    os.chdir(args.directory)

    for directory in os.listdir(args.directory):
        if os.path.isdir(directory) and re.fullmatch("[a-zA-Z0-9]{4}_[0-9]", directory):  # is pdb folder
            make_cross_links(os.path.join(args.directory, directory))

