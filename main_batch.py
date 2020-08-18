
import sys
import os
import subprocess

PWD = "/cs/labs/dina/tomer.cohen13/NR_H_Protein_Martin_1000/"
INTRO = "#!/bin/tcsh\n" \
        "#SBATCH --mem=4000m\n" \
        "#SBATCH -c1\n" \
        "#SBATCH --time=24:0:0\n" \
        "cd /cs/labs/dina/tomer.cohen13\n" \
        "python3 /cs/usr/tomer.cohen13/Desktop/nanobodies/nanobodies_script.py " + PWD


if __name__ == '__main__':
    for file in os.listdir(sys.argv[1]):
        if file.endswith('.pdb'):
            script_name = PWD + file + ".sh"
            with open(script_name, 'w') as f:
                line = INTRO + file.split(".")[0] + " " + file
                f.write(line)
            subprocess.run("sbatch " + script_name, shell=True)
