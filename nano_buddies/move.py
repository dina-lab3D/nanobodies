import os
import subprocess
import pandas as pd
from tqdm import tqdm
import re
import pickle

if __name__ == '__main__':

    dir = "/cs/usr/tomer.cohen13/lab/NN/RosettaFasta"
    os.chdir(dir)

    for pdb_dir in os.listdir(os.getcwd()):
        os.chdir(pdb_dir)

        with open(pdb_dir+'.fa', 'r') as file:
            # read a list of lines into data
            data = file.readlines()

        data[0] = '>heavy\n'

        # and write everything back
        with open(pdb_dir+'.fa', 'w') as file:
            file.writelines(data)
        os.chdir("..")
    exit()

    dir = "/cs/usr/tomer.cohen13/lab/NN/TestPDBs_backup"
    os.chdir(dir)
    if not os.path.exists("/cs/usr/tomer.cohen13/lab/NN/RosettaFasta"):
        os.mkdir("/cs/usr/tomer.cohen13/lab/NN/RosettaFasta")
    for pdb_dir in os.listdir(os.getcwd()):
        pdb_name = os.path.basename(pdb_dir)
        if not os.path.exists("/cs/labs/dina/tomer.cohen13/NN/RosettaFasta/" + pdb_name):
            os.mkdir("/cs/labs/dina/tomer.cohen13/NN/RosettaFasta/" + pdb_name)
        os.chdir(pdb_dir)
        subprocess.run("cp " + pdb_name + ".fa /cs/labs/dina/tomer.cohen13/NN/RosettaFasta/" + pdb_name + "/" + pdb_name + ".fa", shell=True)
        subprocess.run("cp " + pdb_name + ".pdb /cs/labs/dina/tomer.cohen13/NN/RosettaFasta/" + pdb_name + "/" + pdb_name + ".pdb", shell=True)
        subprocess.run("cp ref.pdb /cs/labs/dina/tomer.cohen13/NN/RosettaFasta/" + pdb_name + "/ref.pdb", shell=True)
        subprocess.run("cp ref_cdr1.pdb /cs/labs/dina/tomer.cohen13/NN/RosettaFasta/" + pdb_name + "/ref_cdr1.pdb", shell=True)
        subprocess.run("cp ref_cdr2.pdb /cs/labs/dina/tomer.cohen13/NN/RosettaFasta/" + pdb_name + "/ref_cdr2.pdb", shell=True)
        subprocess.run("cp ref_cdr3.pdb /cs/labs/dina/tomer.cohen13/NN/RosettaFasta/" + pdb_name + "/ref_cdr3.pdb", shell=True)
        os.chdir("..")

    exit()


    if not os.path.exists("/cs/labs/dina/tomer.cohen13/NN/TestPDBs"):
        os.mkdir("/cs/labs/dina/tomer.cohen13/NN/TestPDBs")

    names_file = "/cs/usr/tomer.cohen13/lab/NN/NanoNetPDBs/test_names_save_2.pkl"
    with open(names_file, "rb") as input_file:
        names_list = pickle.load(input_file)

    os.chdir('/cs/usr/tomer.cohen13/lab/NN/NanoNetPDBs')

    for pdb_dir in names_list:

        pdb_name = os.path.basename(pdb_dir)
        if not os.path.exists("/cs/labs/dina/tomer.cohen13/NN/TestPDBs/" + pdb_name):
            os.mkdir("/cs/labs/dina/tomer.cohen13/NN/TestPDBs/" + pdb_name)
        os.chdir(pdb_dir)
        subprocess.run("cp " + pdb_name + ".fa /cs/labs/dina/tomer.cohen13/NN/TestPDBs/" + pdb_name + "/" + pdb_name + ".fa", shell=True)
        subprocess.run("cp " + pdb_name + ".pdb /cs/labs/dina/tomer.cohen13/NN/TestPDBs/" + pdb_name + "/" + pdb_name + ".pdb", shell=True)

        for i in range(5):
            subprocess.run("cp model_" + str(i) + ".pdb /cs/labs/dina/tomer.cohen13/NN/TestPDBs/" + pdb_name + "/model_" + str(i) + ".pdb", shell=True)
            subprocess.run("cp loop_" + str(i) + ".pdb /cs/labs/dina/tomer.cohen13/NN/TestPDBs/" + pdb_name + "/loop_" + str(i) + ".pdb", shell=True)

        subprocess.run("cp ref.pdb /cs/labs/dina/tomer.cohen13/NN/TestPDBs/" + pdb_name + "/ref.pdb", shell=True)
        subprocess.run("cp ref_cdr1.pdb /cs/labs/dina/tomer.cohen13/NN/TestPDBs/" + pdb_name + "/ref_cdr1.pdb", shell=True)
        subprocess.run("cp ref_cdr2.pdb /cs/labs/dina/tomer.cohen13/NN/TestPDBs/" + pdb_name + "/ref_cdr2.pdb", shell=True)
        subprocess.run("cp ref_cdr3.pdb /cs/labs/dina/tomer.cohen13/NN/TestPDBs/" + pdb_name + "/ref_cdr3.pdb", shell=True)

        subprocess.run("cp scores.txt /cs/labs/dina/tomer.cohen13/NN/TestPDBs/" + pdb_name + "/scores.txt", shell=True)

        os.chdir("../..")

    exit()


    os.chdir("/cs/usr/tomer.cohen13/lab/NN/NanoNetPDBs")
    for folder_n in os.listdir(os.getcwd()):
        if os.path.isdir(folder_n) and re.fullmatch('[0-9]+', folder_n):
            os.chdir(folder_n)
            for pdb in tqdm(os.listdir(os.getcwd())):
                os.chdir(pdb)
                if os.path.exists(pdb):
                    os.chdir(pdb)
                    subprocess.run("mv * ..", shell=True)
                    os.chdir("..")
                    subprocess.run("rm -rf " + pdb, shell=True)
                os.chdir("..")
            os.chdir("..")

            #     if os.path.exists("nanobodies.tar"):
            #         if os.path.exists("model_0.pdb"):
            #             df = pd.read_csv("top_models_rmsd.csv", header=None)[1]
            #             for i in range(10):
            #                 os.rename("model_" + str(i) + ".pdb", df[i])
            #         subprocess.run("tar -xf nanobodies.tar", shell=True)
            #         subprocess.run("rm -f nanobodies.tar", shell=True)
            #     os.chdir("..")
            # os.chdir("..")

    # faild_pdbs = pd.read_csv("/cs/usr/tomer.cohen13/lab/NN/NanoNetPDBs/nn_input_failed_pdbs.csv")
    # os.mkdir("/cs/usr/tomer.cohen13/lab/NN/NanoNetPDBs/failed")
    # for pdb, folder in zip(faild_pdbs["PDB"], faild_pdbs["FOLDER"]):
    #     subprocess.run("cp /cs/usr/tomer.cohen13/lab/NN/NanoNetPDBs/" + str(folder) + "/" + pdb + "/" + pdb + ".pdb /cs/usr/tomer.cohen13/lab/NN/NanoNetPDBs/failed/" + pdb + ".pdb", shell=True)
