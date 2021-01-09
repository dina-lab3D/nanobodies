
import argparse
import os
import subprocess
import re



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="directory path containing the pdb directories")
    args = parser.parse_args()
    os.chdir(args.directory)

    for directory in os.listdir(args.directory):
        #  if the folder is pdb folder
        if os.path.isdir(directory) and re.fullmatch("[a-zA-Z0-9]{4}_[0-9]", directory):
            os.chdir(os.path.join(os.getcwd(), directory))
            if os.path.exists( "pydock_results"):
                os.chdir("..")
                continue
            os.mkdir("pydock_results")
            for file in os.listdir(os.getcwd()):
                if (re.fullmatch("docking_(?:loop|model)_[0-9].res", file) or
                re.fullmatch("params_(?:loop|model)_[0-9].txt", file) or
                re.fullmatch("soap_score_(?:loop|model)_[0-9].res", file) or
                re.fullmatch("trans_(?:loop|model)_[0-9]", file) or
                    (file =="best.res" or file =="cat_soap_scores.res" or file =="cdrs3" or file =="soap_score_cluster.res" or
                     file =="ref_scores.csv" or file =="patch_dock.log" or file =="log" or file =="frame" or file =="dock_data.csv" or
                     file =="dist_constraints.read" or file =="dist_constraints")):

                    subprocess.run("mv " + os.path.join(os.getcwd(), file) + " " + os.path.join(os.getcwd(), "pydock_results"),shell=True)
            os.chdir("..")



