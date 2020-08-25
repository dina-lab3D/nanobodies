
import sys
import os
import subprocess

# max memory for each batch
MEMORY = "8000m"

# max time for each batch
TIME = "48:0:0"

# the begining of the script for cluster
INTRO = "#!/bin/tcsh\n" \
        "#SBATCH --mem="+ MEMORY +"\n" \
        "#SBATCH -c1\n" \
        "#SBATCH --time=" + TIME + "\n" \
        "cd /cs/labs/dina/tomer.cohen13\n" \
        "python3 /cs/usr/tomer.cohen13/Desktop/nanobodies/nanobodies_script.py "


if __name__ == '__main__':
    # make sure that you go to cluster (like hm) before calling this program!

    pwd = sys.argv[1]  # working directory
    for file in os.listdir(pwd):
        if file.endswith('.pdb'):  # goes over all pdb files in that directory
            script_name = pwd + "/" + file + ".sh"  # script file
            with open(script_name, 'w') as f:
                line = INTRO + pwd + "/" + file.split(".")[0] + " " + file 
                f.write(line)
            subprocess.run("sbatch " + script_name, shell=True)  # sends script to the cluster
