
import sys
import os
import subprocess


MODE = 'loop_restraints'  # "loop_restraints" or "model_nanobody

# max memory for each batch
MEMORY = "7000m"

# max time for each batch
TIME = "24:0:0"

# NanoModelScript.py path
SCRIPT_PATH = "/cs/labs/dina/tomer.cohen13/nanobodies/nano_buddies/NanoModelScript.py"

# -t for testing mode, -l <int> for number of loops to model
FLAGS = "-t -l 100"

# LoopRestraints.py path
RESTRAINTS_PATH = "/cs/labs/dina/tomer.cohen13/nanobodies/nano_buddies/LoopRestraints.py"

# trained nano_net model path
NANO_NET_MODEL_PATH = "/cs/labs/dina/tomer.cohen13/NN/NanoNetPDBs/save_2"


# the begining of the script for cluster
INTRO = "#!/bin/tcsh\n" \
        "#SBATCH --mem="+ MEMORY +"\n" \
        "#SBATCH -c1\n" \
        "#SBATCH --time=" + TIME + "\n"


if __name__ == '__main__':
    # make sure that you go to cluster (like hm) before calling this program!

    pwd = sys.argv[1]  # working directory
    os.chdir(pwd)


    # for moidelNanobody script
    if MODE == "model_nanobody":
        for file in os.listdir(os.getcwd()):
            if file.endswith('.pdb'):  # goes over all pdb files in that directory
                folder_name = file.split(".")[0]  # the new pdb folder
                if not os.path.isdir(folder_name):
                    os.mkdir(folder_name)

                os.chdir(folder_name)
                script_name = folder_name + ".sh"  # script file

                with open(script_name, 'w') as f:
                    f.write(INTRO)
                    f.write("cd " + os.getcwd() + "\n")

                    # get chain H (antibody)
                    f.write("~dina/utils/getChain.Linux H " + os.getcwd() + ".pdb" + " > " + os.getcwd() + "/ref.pdb")

                    # fasta script
                    f.write("~dina/utils/pdb2fasta " + os.getcwd() + "/ref.pdb" + " > " + os.getcwd() + "/" + folder_name + ".fa")

                    # run the script that creates the loops models
                    f.write("~dina/modeller9.18/bin/modpy.sh python3 /cs/labs/dina/tomer.cohen13/nanobodies/scripts/modelNanobody.py " + FLAGS + " " + os.getcwd() +
                            "/" + folder_name + ".fa")

                    # move the pdb file
                    f.write("mv " + os.getcwd() + ".pdb " + os.getcwd())
                subprocess.run("sbatch " + script_name, shell=True)  # sends script to the cluster

                os.chdir("..")


    # for LoopRestraints script
    if MODE == "loop_restraints":
        for pdb_dir in os.listdir(os.getcwd()):
            os.chdir(pdb_dir)
            script_name = pdb_dir + ".sh"  # script file
            with open(script_name, 'w') as f:
                f.write(INTRO)
                f.write("cd " + os.getcwd() + "\n")
                f.write("module load tensorflow/2.0.0" + "\n")
                f.write("source  ~/venv/plot/bin/activate.csh" + "\n")
                f.write("~dina/modeller9.18/bin/modpy.sh python3 " + RESTRAINTS_PATH + " " + os.getcwd() + " " + NANO_NET_MODEL_PATH)
            subprocess.run("sbatch " + script_name, shell=True)  # sends script to the cluster
            os.chdir("..")

    # for rosetta modeling (no loop modeling)
    if MODE == 'rosetta_model':
        for pdb_dir in os.listdir(os.getcwd()):
            os.chdir(pdb_dir)
            script_name = pdb_dir + ".sh"  # script file
            with open(script_name, 'w') as f:
                f.write(INTRO)
                f.write("cd " + os.getcwd() + "\n")
                f.write("setenv ROSETTA /cs/labs/dina/tomer.cohen13/Rosetta3.12\n")
                f.write("setenv ROSETTA3_DB $ROSETTA/main/database\n")
                f.write("setenv ROSETTA_BIN $ROSETTA/main/source/bin\n")
                f.write("setenv PATH $PATH':'$ROSETTA_BIN\n")
                f.write("antibody.linuxgccrelease -exclude_homologs true -vhh_only -out:file:scorefile scores.txt -fasta " + pdb_dir+".fa | tee grafting.log\n")
                f.write("cd grafting\n")
                f.write("rm -f debug*\n")
                f.write("rm -f orientation*\n")
                f.write("rm -f frh*.pdb*\n")
                f.write("rm -f frl*\n")
            subprocess.run("sbatch " + script_name, shell=True)  # sends script to the cluster
            os.chdir("..")




