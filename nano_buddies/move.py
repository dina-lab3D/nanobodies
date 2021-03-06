import os
import subprocess
import pandas as pd
from tqdm import tqdm
import re
import pickle

if __name__ == '__main__':








    # dir = "/cs/usr/tomer.cohen13/lab/NN/RosettaFasta"
    # # flags_path = "/cs/labs/dina/tomer.cohen13/NN/abH3.flags"
    # os.chdir(dir)
    # for pdb in os.listdir(os.getcwd()):
    #
    #     os.chdir(pdb)
    #     subprocess.run("rm -f slurm-*", shell=True)
    #     subprocess.run("rm -f *.sh", shell=True)
    #     os.chdir("..")
    #     continue
    #
    #     if pdb == "1YC7_1" or pdb == "summery_rosetta_nn" or pdb == "plots_rosetta_nn":
    #         continue
    #     os.chdir(pdb)
    #     for i in range(1, 201):
    #         folder = "H3_modeling"
    #         code = "model-0.relaxed_%04d.pdb" % i
    #         model_path = os.path.join(folder, code)
    #         if not os.path.exists(model_path):
    #             print("failed!!!!!")
    #             print(pdb)
    #             print(i)
    #     print(pdb)
    #     os.chdir("..")
    # exit(0)

    # dir = "/cs/usr/tomer.cohen13/lab/NN/TestPDBs_backup"
    # os.chdir(dir)
    # if not os.path.exists("/cs/usr/tomer.cohen13/lab/NN/RosettaFasta"):
    #     os.mkdir("/cs/usr/tomer.cohen13/lab/NN/RosettaFasta")
    # for pdb_dir in os.listdir(os.getcwd()):
    #     pdb_name = os.path.basename(pdb_dir)
    #     if not os.path.exists("/cs/labs/dina/tomer.cohen13/NN/RosettaFasta/" + pdb_name):
    #         os.mkdir("/cs/labs/dina/tomer.cohen13/NN/RosettaFasta/" + pdb_name)
    #     os.chdir(pdb_dir)
    #     subprocess.run("cp " + pdb_name + ".fa /cs/labs/dina/tomer.cohen13/NN/RosettaFasta/" + pdb_name + "/" + pdb_name + ".fa", shell=True)
    #     subprocess.run("cp " + pdb_name + ".pdb /cs/labs/dina/tomer.cohen13/NN/RosettaFasta/" + pdb_name + "/" + pdb_name + ".pdb", shell=True)
    #     subprocess.run("cp ref.pdb /cs/labs/dina/tomer.cohen13/NN/RosettaFasta/" + pdb_name + "/ref.pdb", shell=True)
    #     subprocess.run("cp ref_cdr1.pdb /cs/labs/dina/tomer.cohen13/NN/RosettaFasta/" + pdb_name + "/ref_cdr1.pdb", shell=True)
    #     subprocess.run("cp ref_cdr2.pdb /cs/labs/dina/tomer.cohen13/NN/RosettaFasta/" + pdb_name + "/ref_cdr2.pdb", shell=True)
    #     subprocess.run("cp ref_cdr3.pdb /cs/labs/dina/tomer.cohen13/NN/RosettaFasta/" + pdb_name + "/ref_cdr3.pdb", shell=True)
    #     os.chdir("..")
    #


    if not os.path.exists("/cs/labs/dina/tomer.cohen13/NN/ModellerFasta"):
        os.mkdir("/cs/labs/dina/tomer.cohen13/NN/ModellerFasta")

    names_file = "/cs/labs/dina/tomer.cohen13/NN/NanoNetPDBs/NanoNet_arrays/test_names_NanoNet_model.pkl"
    with open(names_file, "rb") as input_file:
        names_list = pickle.load(input_file)

    os.chdir('/cs/usr/tomer.cohen13/lab/NN/NanoNetPDBs')

    for pdb_dir in names_list:

        pdb_name = os.path.basename(pdb_dir)
        if not os.path.exists("/cs/labs/dina/tomer.cohen13/NN/ModellerFasta/" + pdb_name):
            os.mkdir("/cs/labs/dina/tomer.cohen13/NN/ModellerFasta/" + pdb_name)
        os.chdir(pdb_dir)
        subprocess.run("cp " + pdb_name + ".fa /cs/labs/dina/tomer.cohen13/NN/ModellerFasta/" + pdb_name + "/" + pdb_name + ".fa", shell=True)
        subprocess.run("cp " + pdb_name + ".pdb /cs/labs/dina/tomer.cohen13/NN/ModellerFasta/" + pdb_name + "/" + pdb_name + ".pdb", shell=True)

        # for i in range(5):
        subprocess.run("cp model_" + str(0) + ".pdb /cs/labs/dina/tomer.cohen13/NN/ModellerFasta/" + pdb_name + "/model_" + str(0) + ".pdb", shell=True)


        subprocess.run("cp ref.pdb /cs/labs/dina/tomer.cohen13/NN/ModellerFasta/" + pdb_name + "/ref.pdb", shell=True)
        subprocess.run("cp ref_cdr1.pdb /cs/labs/dina/tomer.cohen13/NN/ModellerFasta/" + pdb_name + "/ref_cdr1.pdb", shell=True)
        subprocess.run("cp ref_cdr2.pdb /cs/labs/dina/tomer.cohen13/NN/ModellerFasta/" + pdb_name + "/ref_cdr2.pdb", shell=True)
        subprocess.run("cp ref_cdr3.pdb /cs/labs/dina/tomer.cohen13/NN/ModellerFasta/" + pdb_name + "/ref_cdr3.pdb", shell=True)

        # subprocess.run("cp scores.txt /cs/labs/dina/tomer.cohen13/NN/RosettaFasta/" + pdb_name + "/scores.txt", shell=True)

        os.chdir("../..")


    exit()

    dir = "/cs/usr/tomer.cohen13/lab/NN/RosettaFasta"
    # flags_path = "/cs/labs/dina/tomer.cohen13/NN/abH3.flags"
    os.chdir(dir)

    for pdb_dir in os.listdir(os.getcwd()):
        if pdb_dir != "1YC7_1" and pdb_dir != "6QN8_1":
            os.chdir(pdb_dir)
            # subprocess.run("rm -f ROSETTA_CRASH.log" ,shell=True)
            # subprocess.run("rm -f h3_nanonet_modeling-0.log", shell=True)
            # subprocess.run("rm -f slurm*", shell=True)
            # # subprocess.run("rm -f abH3.flags", shell=True)
            # subprocess.run("rm -rf H3_NanoNet_modeling", shell=True)
            # os.mkdir("H3_modeling")
            subprocess.run("rm -f H3_NanoNet_modeling_scores.fasc", shell=True)
            # subprocess.run("cp /cs/labs/dina/tomer.cohen13/NN/abH3.flags /cs/labs/dina/tomer.cohen13/NN/RosettaFasta/" + pdb_dir,shell=True)

        # with open(pdb_dir+'.fa', 'r') as file:
        #     # read a list of lines into data
        #     data = file.readlines()
        #
        # data[0] = '>heavy\n'
        #
        # # and write everything back
        # with open(pdb_dir+'.fa', 'w') as file:
        #     file.writelines(data)
        #
        #
            os.chdir("..")
    exit()


    # os.chdir("/cs/usr/tomer.cohen13/lab/NN/NanoNetPDBs")
    # for folder_n in os.listdir(os.getcwd()):
    #     if os.path.isdir(folder_n) and re.fullmatch('[0-9]+', folder_n):
    #         os.chdir(folder_n)
    #         print(folder_n)
    #         for pdb in tqdm(os.listdir(os.getcwd())):
    #             os.chdir(pdb)
    #             # if os.path.exists(pdb):
    #             #     os.chdir(pdb)
    #             #     subprocess.run("mv * ..", shell=True)
    #             #     os.chdir("..")
    #             #     subprocess.run("rm -rf " + pdb, shell=True)
    #             subprocess.run("rm -f slurm-*", shell=True)
    #             os.chdir("..")
    #         os.chdir("..")

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


    # os.chdir("/cs/usr/tomer.cohen13/lab/nanobodies/COVID_19/nanobodies")
    # for folder_n in os.listdir(os.getcwd()):
    #     if os.path.isdir(folder_n):
    #         os.chdir(folder_n)
    #         subprocess.run("rm -f nanobody_trans*", shell=True)
    #         os.chdir("..")