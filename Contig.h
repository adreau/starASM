#ifndef DEF_CONTIG
#define DEF_CONTIG

#include <string>
#include <vector>

#include "Globals.h"
#include "Interval.h"
#include "Molecule.h"


// Input contigs are usually are split by splitASM
// This is part of a split contig
struct ContigPart: public Interval {

  std::vector < std::string > barcodes_beg;
  std::vector < std::string > barcodes_end;

  ContigPart (unsigned long int pos_beg, unsigned long int pos_end): Interval(pos_beg, pos_end) {}

  void add_beg_molecule (std::string &barcode) {
    barcodes_beg.push_back(barcode);
  }

  void add_end_molecule (std::string &barcode) {
    barcodes_end.push_back(barcode);
  }

  void sort_barcodes () {
    std::sort(barcodes_beg.begin(), barcodes_beg.end());
    std::sort(barcodes_end.begin(), barcodes_end.end());
  }

  void isNeighbourSize (ContigPart &ctgPart, std::vector < int > &arcs);
};


struct Contig {

  std::string name;
  std::vector < ContigPart > contigParts;

  void addPart (unsigned long int pos_beg, unsigned long int pos_end) {
    contigParts.emplace_back(pos_beg, pos_end);
  }

  Contig (std::string &n, unsigned long int pos_beg, unsigned long int pos_end): name(n) {
    addPart(pos_beg, pos_end);
  }

  void sort_barcodes () {
    for (ContigPart &contigPart: contigParts) {
      contigPart.sort_barcodes();
    }
  }
};

#endif
