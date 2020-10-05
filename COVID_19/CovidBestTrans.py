import argparse
import os
import pandas as pd
import numpy as np

NUM_TRANS = 100

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="directory path containing the pdb dock data")
    args = parser.parse_args()

    data_df = pd.read_csv(os.path.join(args.directory, "dock_data.csv"))
    data_df["soap_score"] = np.asfarray(data_df["soap_score"], float)  # cast to float
    top_trans = pd.DataFrame.sort_values(data_df, by="soap_score")[0:NUM_TRANS]["trans"]
    top_trans.to_csv(os.path.join(args.directory, "best_" + str(NUM_TRANS) + "_trans"), header=False, index=False)
