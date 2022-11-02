#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <iostream>
#include "Contig.h"

using namespace std;

Contig::Contig(string &ctg_split_name, string &ctg_orig, int ctg_orig_pos_beg, int ctg_orig_pos_end){

    name = ctg_split_name;
    origin = ctg_orig;
    pos_beg = ctg_orig_pos_beg;
    pos_end = ctg_orig_pos_end;

}

void Contig::add_beg_molecule(Molecule &mol){

    molecules_beg.push_back(mol);
    barcodes_beg.push_back(mol.barcode);
   
}

void Contig::add_mid_molecule(Molecule &mol){

    molecules_mid.push_back(mol);
    barcodes_mid.push_back(mol.barcode);
    
}

void Contig::add_end_molecule(Molecule &mol){

    molecules_end.push_back(mol);
    barcodes_end.push_back(mol.barcode);
    
}

void Contig::sort_barcodes(){

    sort(barcodes_beg.begin(), barcodes_beg.end());

    if(pos_end-pos_beg>200000)
        sort(barcodes_mid.begin(), barcodes_mid.end());

    sort(barcodes_end.begin(), barcodes_end.end());

}

int intersectMoleculesSize(vector<string> &v1, vector<string> &v2, int condition){


    vector<string> commun_molecules;

    set_intersection(v1.begin(),v1.end(),
                     v2.begin(),v2.end(),
                     back_inserter(commun_molecules));
    
    switch (condition)
    {
    case 1:
        if ((commun_molecules.size()>0) && (commun_molecules.size() >= v1.size()*0.8) && (commun_molecules.size() >= v2.size()*0.8))
            return commun_molecules.size();
        else return 0;
        break;
    case 2: 
        if ((commun_molecules.size()>0) && ((commun_molecules.size() >= v1.size()*0.8) || (commun_molecules.size() >= v2.size()*0.8)) )
            return commun_molecules.size();
        else return 0;
        break;
    case 3: 
        if ((commun_molecules.size()>0) && ((commun_molecules.size() >= v1.size()*0.6) && (commun_molecules.size() >= v2.size()*0.6)) )
            return commun_molecules.size();
        else return 0;
        break;
    case 4: 
        if ((commun_molecules.size()>0) && ((commun_molecules.size() >= v1.size()*0.6) || (commun_molecules.size() >= v2.size()*0.6)) )
            return commun_molecules.size();
        else return 0;
        break;
    case 5: 
        if ((commun_molecules.size()>0) && ((commun_molecules.size() >= v1.size()*0.4) && (commun_molecules.size() >= v2.size()*0.4)) )
            return commun_molecules.size();
        else return 0;
        break;
    case 6: 
        if ((commun_molecules.size()>0) && ((commun_molecules.size() >= v1.size()*0.4) || (commun_molecules.size() >= v2.size()*0.4)) )
            return commun_molecules.size();
        else return 0;
        break;
    case 7: 
        if ((commun_molecules.size()>0) && ((commun_molecules.size() >= v1.size()*0.2) && (commun_molecules.size() >= v2.size()*0.2)) )
            return commun_molecules.size();
        else return 0;
        break;
    case 8: 
        if ((commun_molecules.size()>0) && ((commun_molecules.size() >= v1.size()*0.2) || (commun_molecules.size() >= v2.size()*0.2)) )
            return commun_molecules.size();
        else return 0;
        break;
    }

    return 0;
}

void Contig::isNeighbourSize(Contig &ctg, int condition, vector<int> &arcs){

    arcs = vector < int > (4);
    // values: "bb", "be", "eb", "ee"

    arcs[0] = intersectMoleculesSize(barcodes_beg,ctg.barcodes_beg, condition);
    arcs[1] = intersectMoleculesSize(barcodes_beg,ctg.barcodes_end, condition);
    arcs[2] = intersectMoleculesSize(barcodes_end,ctg.barcodes_beg, condition);
    arcs[3] = intersectMoleculesSize(barcodes_end,ctg.barcodes_end, condition);

    if ((ctg.origin == "ctg77" && origin == "ctg107") ||
	(ctg.origin == "ctg107" && origin == "ctg77")){
	
        for(int i=0; i<4;i++){
            cout << "test arc:"<< name <<"\t"<< ctg.name<<"\t"<<arcs[i]<<endl;
        }	
    }
}
