#include <string>

#include "Molecule.h"

using namespace std;

int Molecule::getBeginPos(){

    return begin;

}

int Molecule::getEndPos(){

    return end;

}
string &Molecule::getBarcode(){

    return barcode;

}

int Molecule::getReadsNumber(){

    return noReads;

}
