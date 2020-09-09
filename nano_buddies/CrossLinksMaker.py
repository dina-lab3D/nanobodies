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

XLM_TOOLS = "/cs/labs/dina/tomer.cohen13/XLM-Tools-master/xlmtools.v1.0.py "
INTERFACE_XL = "/cs/staff/dina/libs/imp_build/bin/interface_cross_links "


def make_cross_links(pdf_folder):

    os.chdir(pdf_folder)
    pdb_name = os.path.basename(pdf_folder) + ".pdb"
    with open("xl_script.sh", 'w') as script_file:
        script_file.write(INTRO)
        script_file.write("cd " + os.getcwd() + "\n")
        script_file.write("module load opencv\nsetenv CGAL_DIR /cs/labs/dina/dina/libs/CGAL\n")
        script_file.write(INTERFACE_XL + " antigen.pdb ref.pdb 25 -n -e\n")
        script_file.write("python3 " + XLM_TOOLS + "-xl_list " + os.path.join(pdf_folder, "cxms_all.dat") + " -pdb " + pdb_name + " -jwalk -depth\n")

    subprocess.run("sbatch xl_script.sh", shell=True)
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
        if os.path.isdir(directory) and re.fullmatch("[a-zA-Z0-9]{4}_[0-9]", directory):
            make_cross_links(os.path.join(args.directory, directory))

