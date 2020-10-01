
import sys
import os
import subprocess

# max memory for each batch
MEMORY = "10000m"

# max time for each batch
TIME = "8:0:0"

SCRIPT_PATH = "/cs/labs/dina/tomer.cohen13/nanobodies/nano_buddies/nanobodies_script.py"

# the begining of the script for cluster
INTRO = "#!/bin/tcsh\n" \
        "#SBATCH --mem="+ MEMORY +"\n" \
        "#SBATCH -c1\n" \
        "#SBATCH --time=" + TIME + "\n"


if __name__ == '__main__':
    # make sure that you go to cluster (like hm) before calling this program!

    pwd = sys.argv[1]  # working directory
    for file in os.listdir(pwd):
        if file.endswith('.pdb'):  # goes over all pdb files in that directory
            current_folder_path = pwd + "/" + file.split(".")[0]

            if not os.path.isdir(current_folder_path):
                os.mkdir(current_folder_path)

            os.chdir(current_folder_path)
            script_name = file + ".sh"  # script file
            with open(script_name, 'w') as f:
                #line = INTRO + pwd + "/" + file.split(".")[0] + " " + file
                f.write(INTRO)
                line = "cd " + os.getcwd() + "\n"
                f.write(line)
                line = "python3 " + SCRIPT_PATH + " " + os.getcwd() + " " + file
                f.write(line)
            subprocess.run("sbatch " + script_name, shell=True)  # sends script to the cluster
            os.chdir("../..")


