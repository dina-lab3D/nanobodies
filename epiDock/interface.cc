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
    mol1.readPDBfile(molFile1, PDB::ChainSelector(chains1));
    molFile1.close();
    std::ifstream molFile2(argv[i]);
    mol2.readAllPDBfile(molFile2, PDB::ChainSelector(chains2));
    molFile2.close();
    MolecularInterface molInterface(mol1, mol2, interfaceThr);
    molInterface.getReceptorInterfaceResidues(epitope);
  }

  std::cout << epitope.size() << std::endl;
  std::ofstream ofile("epi.csv");
  int numberOfChains  = 1;// TODO: change back to 1, add as option later
  int range = epitope.size()/numberOfChains;
  int molNum = argc -4;
  for(int i=0; i<range; i++) {
    int counter = epitope[i];// first chain
    for(int chain = 1; chain < numberOfChains; chain++) {
      counter += epitope[i+range*chain];
    }
    float ratio = (float)counter/molNum;
    ofile << ratio << ", ";
  }
  ofile << std::endl;
  ofile.close();
  return 0;
}
