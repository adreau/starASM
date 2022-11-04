#ifndef DEF_MOLECULE
#define DEF_MOLECULE

#include <string>
#include "Interval.h"

struct Molecule: public Interval {
  
  Molecule(unsigned long int b, unsigned long int e, std::string &bc, int nr)
      : Interval(b, e), barcode(bc), noReads(nr) {}

  int noReads;

  std::string barcode;
};


#endif
