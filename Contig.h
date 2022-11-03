#ifndef DEF_CONTIG
#define DEF_CONTIG

#include <string>
#include <vector>

#include "Molecule.h"


class Contig
{

  public:

    std::string name;
    std::string origin;
    int pos_beg,pos_end;
    std::vector<std::string> barcodes_beg;
    std::vector<std::string> barcodes_mid;
    std::vector<std::string> barcodes_end;

    Contig(std::string &name, std::string &origin, int pos_beg, int pos_end);

    void add_beg_molecule(Molecule &mol);
    void add_mid_molecule(Molecule &mol);
    void add_end_molecule(Molecule &mol);
    void sort_barcodes();
    
    void isNeighbourSize(Contig &ctg, int condition, std::vector < int > &arcs);

};

#endif
