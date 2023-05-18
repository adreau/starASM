#ifndef molecule_h
#define molecule_h

#include <iostream>
#include <string>

#include "interval.h"


struct Molecule: public Interval {
  unsigned int chrid;
  std::string barcode;
  unsigned int n_reads;
  Molecule (unsigned int c, unsigned long s, unsigned long e, const std::string &b, unsigned int n): Interval(s, e), chrid(c), barcode(b), n_reads(n) {}
  friend bool operator< (const Molecule& lhs, const Molecule& rhs);
  friend std::ostream& operator<<(std::ostream& os, const Molecule& molecule);
};

using Molecules = std::vector < Molecule >;

#endif
