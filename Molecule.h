#ifndef DEF_MOLECULE
#define DEF_MOLECULE

#include <string>

using namespace std;

class Molecule {
  
 public:

  Molecule(int &b, int &e, string &bc, int &nR)
      : begin(b), end(e), barcode(bc), noReads(nR) {}

  int begin, end, noReads;

  string barcode;

  int getBeginPos();
  int getEndPos();
  string &getBarcode();
  int getReadsNumber();




};


#endif
