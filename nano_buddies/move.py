
import os
import subprocess
import pandas as pd
from tqdm import tqdm

if __name__ == '__main__':
    os.chdir("/cs/usr/tomer.cohen13/lab/NN/NanoNetPDBs")
    for folder_n in os.listdir(os.getcwd()):
        if os.path.isdir(folder_n) and folder_n != 'failed':
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

#                if os.path.exists("nanobodies.tar"):
#
#                    if os.path.exists("model_0.pdb"):
#                        df = pd.read_csv("top_models_rmsd.csv", header=None)[1]
#                        for i in range(10):
#                            os.rename("model_" + str(i) + ".pdb", df[i])
#                    subprocess.run("tar -xf nanobodies.tar", shell=True)
#                    subprocess.run("rm -f nanobodies.tar", shell=True)
#                os.chdir("..")
#            os.chdir("..")









    # faild_pdbs = pd.read_csv("/cs/usr/tomer.cohen13/lab/NN/NanoNetPDBs/nn_input_failed_pdbs.csv")
    # os.mkdir("/cs/usr/tomer.cohen13/lab/NN/NanoNetPDBs/failed")
    # for pdb, folder in zip(faild_pdbs["PDB"], faild_pdbs["FOLDER"]):
    #     subprocess.run("cp /cs/usr/tomer.cohen13/lab/NN/NanoNetPDBs/" + str(folder) + "/" + pdb + "/" + pdb + ".pdb /cs/usr/tomer.cohen13/lab/NN/NanoNetPDBs/failed/" + pdb + ".pdb", shell=True)
