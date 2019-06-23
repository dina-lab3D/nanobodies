#include "MolecularInterface.h"

#include <unistd.h>
#include <stdio.h>
#include <limits.h>

#define CHAIN1_IDX 1
#define CHAIN2_IDX 2
#define INTERFACE_THR_IDX 3

using std::ifstream;

class ChainNoWaterHydrogenSelector : public PDB::Selector {
public:
    // constructor
    ChainNoWaterHydrogenSelector(const std::string &chainID) : chains(chainID) {}

    virtual ~ChainNoWaterHydrogenSelector() {}

    bool operator()(const char *PDBrec) const {
        PDB::WaterHydrogenUnSelector sel;
        if (sel(PDBrec)) {
            PDB::ChainSelector chainSel(chains);
            return chainSel(PDBrec);
        }
        return false;
    }

private:
    std::string chains;
};

int main(int argc, char **argv) {
  if (argc < 5) {
    std::cout << "usage " << argv[0] << " AB H thr pdb1 pdb2 ... " << std::endl;
    return 0;
  }

  std::string chains1 = argv[CHAIN1_IDX];
  std::string chains2 = argv[CHAIN2_IDX];
  float interfaceThr = atof(argv[INTERFACE_THR_IDX]);

  std::vector<double> epitope;

  for(int i=INTERFACE_THR_IDX+1; i<argc; i++) {
    /* read molecules */
    ChemMolecule mol1, mol2;

    std::ifstream molFile1(argv[i]);
    mol1.readAllPDBfile(molFile1, PDB::ChainSelector(chains1));
    molFile1.close();
    std::ifstream molFile2(argv[i]);
    mol2.readAllPDBfile(molFile2, PDB::ChainSelector(chains2));
    molFile2.close();
    MolecularInterface molInterface(mol1, mol2, interfaceThr);
    molInterface.getReceptorInterfaceResidues(epitope);
  }

  std::ofstream ofile("epi.csv");
  int shift = epitope.size()/2;
  for(int i=0; i<shift; i++) {
    ofile << epitope[i]+epitope[i+shift] << ", ";
  }
  ofile << std::endl;
  ofile.close();
  return 0;
}
