#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <iostream>
#include "Contig.h"

using namespace std;

Contig::Contig(string ctg_split_name, string ctg_orig, int ctg_orig_pos_beg, int ctg_orig_pos_end){

    name = ctg_split_name;
    origin = ctg_orig;
    pos_beg = ctg_orig_pos_beg;
    pos_end = ctg_orig_pos_end;

}

void Contig::add_beg_molecule(Molecule mol){

    molecules_beg.push_back(mol);
    barcodes_beg.push_back(mol.barcode);
   
}

void Contig::add_mid_molecule(Molecule mol){

    molecules_mid.push_back(mol);
    barcodes_mid.push_back(mol.barcode);
    
}

void Contig::add_end_molecule(Molecule mol){

    molecules_end.push_back(mol);
    barcodes_end.push_back(mol.barcode);
    
}

void Contig::sort_barcodes(){

    sort(barcodes_beg.begin(), barcodes_beg.end());

    if(pos_end-pos_beg>200000)
        sort(barcodes_mid.begin(), barcodes_mid.end());

    sort(barcodes_end.begin(), barcodes_end.end());

}

int intersectMoleculesSize(vector<string> v1, vector<string> v2, int condition){


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

}

/*vector<int> Contig::isNeighbour(Contig ctg, int condition){

    vector<int> arcs;
    for(int i=0;i<4;i++){
        arcs.push_back(0);
    } // values: "bb", "be", "eb", "ee"

    if(intersectMolecules(barcodes_beg,ctg.barcodes_beg, condition)) arcs[0] = 1;
    if(intersectMolecules(barcodes_beg,ctg.barcodes_end, condition)) arcs[1] = 1;
    if(intersectMolecules(barcodes_end,ctg.barcodes_beg, condition)) arcs[2] = 1;
    if(intersectMolecules(barcodes_end,ctg.barcodes_end, condition)) arcs[3] = 1;


    return arcs;
}*/

vector<int> Contig::isNeighbourSize(Contig ctg, int condition){

    vector<int> arcs;
    for(int i=0;i<4;i++){
        arcs.push_back(0);
    } // values: "bb", "be", "eb", "ee"

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

    /* add mid molecules if short contig => no longer necessary
    if((arcs[0]<5)&&(arcs[1]<5)&&(arcs[2]<5)&&(arcs[3]<5)) {
         
        int length_ctg1 = pos_end - pos_beg;
        int length_ctg2 = ctg.pos_end - ctg.pos_beg;

        if (length_ctg1 < 200000){ //&& (length_ctg2 >= 200000)){

            int mid_beg = intersectMoleculesSize(barcodes_mid,ctg.barcodes_beg);
            int mid_end = intersectMoleculesSize(barcodes_mid,ctg.barcodes_end);
            arcs[0] += mid_beg;
            arcs[1] += mid_end;
            arcs[2] += mid_beg;
            arcs[3] += mid_end;

        }

        if (length_ctg2 < 200000){ //&& (length_ctg1 >= 200000)){

            int beg_mid = intersectMoleculesSize(barcodes_beg,ctg.barcodes_mid);
            int end_mid = intersectMoleculesSize(barcodes_end,ctg.barcodes_mid);
            arcs[0] += beg_mid;
            arcs[1] += beg_mid;
            arcs[2] += end_mid;
            arcs[3] += end_mid;

        }

    } */
        


    return arcs;
}
