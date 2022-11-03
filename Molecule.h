#ifndef DEF_MOLECULE
#define DEF_MOLECULE

#include <string>

class Molecule {
  
 public:

  Molecule(int &b, int &e, std::string &bc, int &nR)
      : begin(b), end(e), barcode(bc), noReads(nR) {}

  int begin, end, noReads;

  std::string barcode;

};


#endif
