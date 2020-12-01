from Bio.PDB import *
import os


def valid_pdb(pdb, test=False):

    if not os.path.exists(os.path.join(pdb, "ref.pdb")):
        return False
    model = PDBParser().get_structure(id=pdb, file=os.path.join(pdb, "ref.pdb"))[0]["H"]
    for i in model.get_residues():
        if not i.has_id("N"):
            print(pdb + ": no side chains")
            return False
        break

    if not test:
        if not os.path.exists(os.path.join(pdb, "model_0.pdb")):
            print(pdb + ": failed modeling")
            return False

    return True


def get_seq(chain):

    aa_residues = []
    seq = ""

    for residue in chain.get_residues():
        aa = residue.get_resname()
        if not is_aa(aa) or not residue.has_id('CA'):
            continue
        elif aa == "UNK":
            seq += "X"
            aa_residues.append(residue)
        else:
            seq += Polypeptide.three_to_one(residue.get_resname())
            aa_residues.append(residue)

    return seq, aa_residues
