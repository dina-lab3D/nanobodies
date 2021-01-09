import sys
import pandas as pd

if __name__ == '__main__':

    file1 = set(pd.read_csv(sys.argv[1], usecols=[0])["PDB"])
    file2 = set(pd.read_csv(sys.argv[2], usecols=[0])["PDB"])

    print(file1 ^ file2)  # symmetric difference

