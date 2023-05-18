#include <cassert>
#include <fstream>
#include <sstream>
#include <iostream>

#include "globals.h"
#include "interval.h"
#include "parsers.h"


void parse_molecule_file (Molecules &molecules) {
  std::ifstream molecule_file (Globals::input_molecules_file_name);
  std::string line, contig, barcode;
  unsigned long start, end;
  unsigned int n_reads;
  if (! molecule_file.is_open()) {
    std::cerr << "Cannot open molecule file.\n";
    exit(EXIT_FAILURE);
  }
  while(getline(molecule_file, line)) {
    std::stringstream ss (line);
    ss >> contig >> start >> end >> barcode >> n_reads;
    assert(Globals::chrids.find(contig) != Globals::chrids.end());
    molecules.emplace_back(Globals::chrids[contig], start, end, barcode, n_reads);
  }
}

void parse_split_file (Contigs &contigs) {
  std::ifstream contig_file(Globals::input_split_file_name);
  std::string contig_line, ctg;
  int pos_beg, pos_end;
  unsigned int n_contig_parts;
  if (! contig_file.is_open()){
      std::cerr << "Error!  Cannot open file '" << Globals::input_split_file_name << "'.\n";
      exit(EXIT_FAILURE);
  }
  contigs.resize(Globals::chrs.size());
  for (n_contig_parts = 0; getline(contig_file, contig_line); ++n_contig_parts) {
    std::stringstream splitstream (contig_line);
    splitstream >> ctg >> pos_beg >> pos_end;
    // BED format is 0-based on the start, and 1-based on the end
    ++pos_beg;
    contigs[Globals::chrids[ctg]].add_part(pos_beg, pos_end);
  }
  std::cerr << n_contig_parts << " contig parts, seen.\n";
}
