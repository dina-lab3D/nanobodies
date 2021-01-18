
import sys
import os
import subprocess

# 'loop_restraints'/'model_nanobody'/'rosetta_model'/'rosetta_loops'
MODE = 'rosetta_loops'

# NanoModelScript.py path
SCRIPT_PATH = "/cs/labs/dina/tomer.cohen13/nanobodies/nano_buddies/NanoModelScript.py"

# LoopRestraints.py path
RESTRAINTS_PATH = "/cs/labs/dina/tomer.cohen13/nanobodies/nano_buddies/LoopRestraints.py"

# trained nano_net model path
NANO_NET_MODEL_PATH = "/cs/labs/dina/tomer.cohen13/NN/NanoNetPDBs/NanoNet_model"

CONST = True

FAILED = ["1YC7_1"]

LONG = ["6QN8_1"]

# the begining of the script for cluster
INTRO = "#!/bin/tcsh\n" \
        "#SBATCH --mem={}\n" \
        "#SBATCH -c1\n" \
        "#SBATCH --time={}\n" #\
        # "#SBATCH --array=1-3\n"


def model_nanobody():
    """

    :return:
    """
    time, memory = "8:0:0", "7000m"
    # -t for testing mode, -l <int> for number of loops to model
    flags = "-t -l 100"
    for file in os.listdir(os.getcwd()):
        if file.endswith('.pdb'):  # goes over all pdb files in that directory
            folder_name = file.split(".")[0]  # the new pdb folder
            if not os.path.isdir(folder_name):
                os.mkdir(folder_name)

            os.chdir(folder_name)
            script_name = folder_name + ".sh"  # script file

            with open(script_name, 'w') as f:
                f.write(INTRO.format(memory, time))
                f.write("cd " + os.getcwd() + "\n")

                # get chain H (antibody)
                f.write("~dina/utils/getChain.Linux H " + os.getcwd() + ".pdb" + " > " + os.getcwd() + "/ref.pdb")

                # fasta script
                f.write("~dina/utils/pdb2fasta " + os.getcwd() + "/ref.pdb" + " > " + os.getcwd() + "/" + folder_name + ".fa")

                # run the script that creates the loops models
                f.write("~dina/modeller9.18/bin/modpy.sh python3 /cs/labs/dina/tomer.cohen13/nanobodies/scripts/modelNanobody.py " + flags + " " + os.getcwd() +"/" + folder_name + ".fa")

                # move the pdb file
                f.write("mv " + os.getcwd() + ".pdb " + os.getcwd())
            subprocess.run("sbatch " + script_name,shell=True)  # sends script to the cluster

            os.chdir("..")


def loop_restraints():
    """

    :return:
    """
    time, memory = "14:0:0", "7000m"

    for pdb_dir in os.listdir(os.getcwd()):
        os.chdir(pdb_dir)
        script_name = pdb_dir + ".sh"  # script file
        with open(script_name, 'w') as f:
            f.write(INTRO.format(memory, time))
            f.write("cd " + os.getcwd() + "\n")
            f.write("module load tensorflow/2.0.0" + "\n")
            f.write("source  ~/venv/plot/bin/activate.csh" + "\n")
            f.write("~dina/modeller9.18/bin/modpy.sh python3 " + RESTRAINTS_PATH + " " + os.getcwd() + " " + NANO_NET_MODEL_PATH)
        subprocess.run("sbatch " + script_name,shell=True)  # sends script to the cluster
        os.chdir("..")


def rosetta_model():
    """

    :return:
    """
    time, memory = "6:0:0", "7000m"

    for pdb_dir in os.listdir(os.getcwd()):
        if pdb_dir not in FAILED:
            continue
            # -detect_disulf false
            # -camelid true
        os.chdir(pdb_dir)
        script_name = pdb_dir + ".sh"  # script file
        with open(script_name, 'w') as f:
            f.write(INTRO.format(memory, time))
            f.write("cd " + os.getcwd() + "\n")
            f.write("setenv ROSETTA /cs/labs/dina/tomer.cohen13/Rosetta\n")
            f.write("setenv ROSETTA3_DB $ROSETTA/main/database\n")
            f.write("setenv ROSETTA_BIN $ROSETTA/main/source/bin\n")
            f.write("setenv PATH $PATH':'$ROSETTA_BIN\n")
            f.write("setenv PATH $PATH':'/cs/labs/dina/tomer.cohen13/Blast/bin\n")
            f.write("antibody.linuxgccrelease -exclude_homologs true -vhh_only -fasta " + pdb_dir + ".fa | tee grafting.log\n")
            f.write("cd grafting\n")
            f.write("rm -f debug*\n")
            f.write("rm -f orientation*\n")
            f.write("rm -f frh*.pdb*\n")
            f.write("rm -f frl*\n")
            f.write("rm -f l*\n")
        subprocess.run("sbatch " + script_name,shell=True)  # sends script to the cluster
        os.chdir("..")


def rosetta_loops():
    """

    :return:
    """
    time, memory, array = "6-0", "4500m", "3"

    for pdb_dir in os.listdir(os.getcwd()):
        if pdb_dir in FAILED or pdb_dir in LONG:
            continue
            # -detect_disulf false
            # -camelid true
        os.chdir(pdb_dir)
        script_name = pdb_dir + ".sh" if not CONST else pdb_dir + "_nanonet.sh" # script file
        with open(script_name, 'w') as f:
            f.write(INTRO.format(memory, time))
            f.write("#SBATCH --array=1-{}\n".format(array))
            f.write("cd " + os.getcwd() + "\n")
            f.write("setenv ROSETTA /cs/labs/dina/tomer.cohen13/Rosetta\n")
            f.write("setenv ROSETTA3_DB $ROSETTA/main/database\n")
            f.write("setenv ROSETTA_BIN $ROSETTA/main/source/bin\n")
            f.write("setenv PATH $PATH':'$ROSETTA_BIN\n")
            f.write("setenv PATH $PATH':'/cs/labs/dina/tomer.cohen13/Blast/bin\n")
            if CONST:
                if not os.path.isdir("H3_NanoNet_modeling"):
                    os.mkdir("H3_NanoNet_modeling")
                f.write("antibody_H3.linuxgccrelease @abH3.flags -s grafting/model-0.relaxed.pdb -nstruct 200 -out:file:scorefile H3_NanoNet_modeling_scores.fasc -out:path:pdb H3_NanoNet_modeling -constraints:cst_file {}_constraints -constraints:cst_weight 1.0 > h3_nanonet_modeling-0.log\n".format(pdb_dir))
            else:
                if not os.path.isdir("H3_modeling"):
                    os.mkdir("H3_modeling")
                f.write("antibody_H3.linuxgccrelease @abH3.flags -s grafting/model-0.relaxed.pdb -nstruct 200 -out:file:scorefile H3_modeling_scores.fasc -out:path:pdb H3_modeling > h3_modeling-0.log\n")
            # for model in range(1,10):
            #     f.write("antibody_H3.linuxgccrelease @abH3.flags -s grafting/model-{}.relaxed.pdb -nstruct 100 > h3_modeling-{}.log\n".format(model, model))
        subprocess.run("sbatch " + script_name,shell=True)  # sends script to the cluster
        os.chdir("..")


if __name__ == '__main__':
    # make sure that you go to cluster (like hm) before calling this program!

    pwd = sys.argv[1]  # working directory
    os.chdir(pwd)

    # for moidelNanobody script
    if MODE == "model_nanobody":
        model_nanobody()

    # for LoopRestraints script
    if MODE == "loop_restraints":
        loop_restraints()

    # for rosetta modeling (no loop modeling)
    if MODE == 'rosetta_model':
        rosetta_model()

    # for rosetta cdr3 loop modeling
    if MODE == 'rosetta_loops':
        rosetta_loops()


