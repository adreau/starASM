#include <iostream>

#include "constants.h"
#include "globals.h"
#include "molecule.h"


bool operator< (const Molecule& lhs, const Molecule& rhs) {
  if (lhs.chrid < rhs.chrid) {
    return true;
  }
  if (lhs.chrid > rhs.chrid) {
    return false;
  }
  return (lhs.start < rhs.start);
}

std::ostream& operator<<(std::ostream& os, const Molecule& molecule) {
  os << Globals::chrs[molecule.chrid] << TAB <<
    molecule.start                    << TAB <<
    molecule.end                      << TAB <<
    molecule.barcode                  << TAB <<
    molecule.n_reads                  << "\n";
  return os;
}