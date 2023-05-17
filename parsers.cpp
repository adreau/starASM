#include <cassert>
#include <vector>
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

void parse_split_file (RefIntervalsSet &refIntervalsSet) {
  std::ifstream joins_file (Globals::input_split_file_name);
  if (! joins_file.is_open()) {
    std::cerr << "Error!  Cannot open file '" << Globals::input_split_file_name << "'" << std::endl;
    exit(EXIT_FAILURE);
  }
  std::string line, part, ref;
  unsigned long int start, end;
  char strand;
  while (getline(joins_file, line)) {
    std::istringstream split_line(line);
    RefIntervals refIntervals;
    while (getline(split_line, part, ';')) {
      if (! part.empty()) {
        std::stringstream split_part(part);
        split_part >> strand >> ref >> start >> end;
        refIntervals.emplace_back(ref, strand, start, end);
      }
    }
    refIntervalsSet.push_back(refIntervals);
  }
}
