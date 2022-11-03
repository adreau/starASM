#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <iostream>
#include "Globals.h"
#include "Contig.h"

Contig::Contig(std::string &ctg_split_name, std::string &ctg_orig, int ctg_orig_pos_beg, int ctg_orig_pos_end){

    name = ctg_split_name;
    origin = ctg_orig;
    pos_beg = ctg_orig_pos_beg;
    pos_end = ctg_orig_pos_end;

}

void Contig::add_beg_molecule(Molecule &mol){

    barcodes_beg.push_back(mol.barcode);
   
}

void Contig::add_end_molecule(Molecule &mol){

    barcodes_end.push_back(mol.barcode);
    
}

void Contig::sort_barcodes(){

    sort(barcodes_beg.begin(), barcodes_beg.end());
    sort(barcodes_end.begin(), barcodes_end.end());

}

int intersectMoleculesSize(std::vector<std::string> &v1, std::vector<std::string> &v2){

    std::vector<std::string> common_molecules;

    set_intersection(v1.begin(),v1.end(),
                     v2.begin(),v2.end(),
                     back_inserter(common_molecules));
    size_t s1 = v1.size();
    size_t s2 = v2.size();
    size_t sc = common_molecules.size();
    if (sc == 0) return 0;
    
    switch (Globals::condition)
    {
    case 1:
        if ((sc >= s1 * 0.8) && (sc >= s2 * 0.8)) return sc;
        return 0;
    case 2: 
        if ((sc >= s1 * 0.8) || (sc >= s2 * 0.8)) return sc;
        return 0;
    case 3: 
        if ((sc >= s1 * 0.6) && (sc >= s2 * 0.6)) return sc;
        return 0;
    case 4: 
        if ((sc >= s1 * 0.6) || (sc >= s2 * 0.6)) return sc;
        return 0;
    case 5: 
        if ((sc >= s1 * 0.4) && (sc >= s2 * 0.4)) return sc;
        return 0;
    case 6: 
        if ((sc >= s1 * 0.4) || (sc >= s2 * 0.4)) return sc;
        return 0;
    case 7: 
        if ((sc >= s1 * 0.2) && (sc >= s2 * 0.2)) return sc;
        return 0;
    case 8: 
        if ((sc >= s1 * 0.2) || (sc >= s2 * 0.2)) return sc;
        return 0;
    }

    return 0;
}

void Contig::isNeighbourSize(Contig &ctg, std::vector<int> &arcs){

    arcs = std::vector < int > (4);
    // values: "bb", "be", "eb", "ee"

    arcs[0] = intersectMoleculesSize(barcodes_beg,ctg.barcodes_beg);
    arcs[1] = intersectMoleculesSize(barcodes_beg,ctg.barcodes_end);
    arcs[2] = intersectMoleculesSize(barcodes_end,ctg.barcodes_beg);
    arcs[3] = intersectMoleculesSize(barcodes_end,ctg.barcodes_end);
}
