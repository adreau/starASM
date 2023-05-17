#ifndef molecule_h
#define molecule_h

#include <iostream>
#include <string>


struct Molecule {
  unsigned int chrid;
  unsigned long start, end;
  std::string barcode;
  unsigned int n_reads;
  Molecule (unsigned int c, unsigned long s, unsigned long e, const std::string &b, unsigned int n): chrid(c), start(s), end(e), barcode(b), n_reads(n) {}
  friend bool operator< (const Molecule& lhs, const Molecule& rhs);
  friend std::ostream& operator<<(std::ostream& os, const Molecule& molecule);
};

using Molecules = std::vector < Molecule >;

#endif
