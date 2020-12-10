import Bio.PDB
import argparse
import pandas as pd


def color(pdb, chains, csv):

    df = pd.read_csv(csv)
    pdb_parser = Bio.PDB.PDBParser()
    model = pdb_parser.get_structure("name", pdb)[0]
    means_df = df.mean(axis=0) * 100
    for chain_id in chains:
        chain = model[chain_id]
        for residue_index, interface in zip(means_df.index, means_df):
            residue = chain[int(residue_index)]
            for atom in residue:
                atom.set_bfactor(interface)
    io = Bio.PDB.PDBIO()
    io.set_structure(model)
    io.save(pdb.replace(".pdb", "_epitops.pdb"))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="pdb to color")
    parser.add_argument("chains", help="chains to color according to csv")
    parser.add_argument("csv", help="csv file containing the interface percentage")
    args = parser.parse_args()

    color(args.pdb, args.chains, args.csv)