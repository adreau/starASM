#ifndef CONTIG_H
#define CONTIG_H

#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>

#include "globals.h"
#include "node.h"
#include "interval.h"
#include "molecule.h"


// Input contigs are usually are split by splitASM
// This is part of a split contig
struct ContigPart: public Interval {

  std::vector < unsigned long int > barcodes_beg;
  std::vector < unsigned long int > barcodes_end;
  std::vector < unsigned long int > all_barcodes;

  ContigPart (unsigned long int pos_beg, unsigned long int pos_end): Interval(pos_beg, pos_end) {}

  void add_beg_molecule (unsigned long int barcode) {
    barcodes_beg.push_back(barcode);
    all_barcodes.push_back(barcode);
  }

  void add_end_molecule (unsigned long int barcode) {
    barcodes_end.push_back(barcode);
    all_barcodes.push_back(barcode);
  }

  void add_other_molecule (unsigned long int barcode) {
    all_barcodes.push_back(barcode);
  }

  void sort_barcodes () {
    std::sort(barcodes_beg.begin(), barcodes_beg.end());
    std::sort(barcodes_end.begin(), barcodes_end.end());
    std::sort(all_barcodes.begin(), all_barcodes.end());
  }
};


struct Contig {

  std::vector < ContigPart > contigParts;

  bool empty () const {
    return contigParts.empty();
  }

  void add_part (unsigned long int pos_beg, unsigned long int pos_end) {
    contigParts.emplace_back(pos_beg, pos_end);
  }

  void sort_barcodes () {
    for (ContigPart &contigPart: contigParts) {
      contigPart.sort_barcodes();
    }
  }
};


using Contigs = std::vector < Contig >;

void add_molecules_to_contigs_extremites (Molecules &molecules, Contigs &contigs);
void intersectMoleculesSize(std::vector < unsigned long int > &b1, std::vector < unsigned long int > &b2, unsigned int &n_inter, double &jaccard);
Contig &get_contig (Contigs &contigs, NodeId &nodeId);
ContigPart &get_contig_part (Contigs &contigs, NodeId &nodeId);


#endif
