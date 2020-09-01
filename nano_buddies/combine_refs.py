import os
import argparse
import subprocess
import pandas as pd

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="directory path to the pdb folder")
    args = parser.parse_args()
    os.chdir(args.directory)

    names = []
    scores = []

    for file in os.listdir(os.getcwd()):
        if file.startswith("no_trans") or file == "soap_score_ref.res":
            soap_score = subprocess.run("grep \"|\" " + file + " | cut -d '|' -f2", shell=True, capture_output=True, universal_newlines=True).stdout
            soap_score = soap_score.replace(" ","").split("\n")[1:-1]  # split to array
            names.append(file.split(".")[0].split("score_")[-1])
            scores.append(soap_score[0])
            os.remove(file)
    df = pd.DataFrame({"name": names, "soap_score": scores})
    df.to_csv("ref_scores.csv", index=False, header=False)
    os.chdir("..")