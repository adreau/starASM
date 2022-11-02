#include <string>

#include "Molecule.h"

int Molecule::getBeginPos(){

    return begin;

}

int Molecule::getEndPos(){

    return end;

}
std::string &Molecule::getBarcode(){

    return barcode;

}

int Molecule::getReadsNumber(){

    return noReads;

}
