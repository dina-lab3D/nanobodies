


import argparse
import pandas as pd
import numpy as np

MAX_N_TEMPLATES = 10

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("csv_file1", help="first summery csv")
    parser.add_argument("csv_file2", help="first summery csv")
    args = parser.parse_args()

    df1 = pd.read_csv(args.csv_file1, usecols=[0,2,5])
    df2 = pd.read_csv(args.csv_file2, usecols=[0,2,5])
	
    first_pdb = True
    with open("compare_rmsd.csv", "w") as sum_file:
        for i in range(len(df1["PDB"])):
            pdb_name = df1["PDB"][i]
            rmsds_df2 = df2[df2["PDB"] == pdb_name]
            if rmsds_df2.empty:
                continue
            print(rmsds_df2)
            df_to_write = pd.DataFrame({"PDB":[pdb_name], "file1_rmsd": [df1["RMSD_BEST_10"][i]], "file2_rmsd": rmsds_df2["RMSD_BEST_10"], "diff_rmsd":[df1["RMSD_BEST_10"][i]-rmsds_df2.iloc[0,1]] , "file1_min_rmsd": [df1["MIN_RMSD"][i]] , "file2_min_rmsd": rmsds_df2["MIN_RMSD"] , "diff_min_rmsd": [df1["MIN_RMSD"][i]- rmsds_df2.iloc[0,2]]})
            df_to_write.to_csv(sum_file, header=first_pdb, index=False)
            first_pdb = False

