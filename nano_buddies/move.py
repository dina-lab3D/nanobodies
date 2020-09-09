
import sys
import os
import subprocess

if __name__ == '__main__':
    os.chdir(sys.argv[1])
    for pdb in os.listdir(sys.argv[1]):
        if pdb.endswith(".pdb"):
            subprocess.run("mv " + pdb + " " + pdb.split(".")[0], shell=True)