#ifndef DEF_CONTIG
#define DEF_CONTIG

#include <string>
#include <vector>

#include "Molecule.h"

using namespace std;


class Contig
{

  public:

    string name;
    string origin;
    int pos_beg,pos_end;
    vector<Molecule> molecules_beg;
    vector<Molecule> molecules_mid;
    vector<Molecule> molecules_end;
    vector<string> barcodes_beg;
    vector<string> barcodes_mid;
    vector<string> barcodes_end;

    Contig(string name, string origin, int pos_beg, int pos_end);

    void add_beg_molecule(Molecule mol);
    void add_mid_molecule(Molecule mol);
    void add_end_molecule(Molecule mol);
    void sort_barcodes();
    
    //vector<int> isNeighbour(Contig ctg, int condition);
    vector<int> isNeighbourSize(Contig ctg, int condition);

};

#endif
