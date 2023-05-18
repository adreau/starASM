#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <map>
#include <stdexcept>    // std::out_of_range

#include "globals.h"
#include "fasta.h"
#include "interval.h"

void extract_subsequence (RefInterval &refInterval, std::string &gap, std::string &output) {
  std::string contig;
  output.clear();
  try {
    contig = Globals::sequences.at(refInterval.ref);
  }
  catch (const std::out_of_range& oor) {
    std::cerr << "Cannot find contig '" << refInterval.ref << "'.\nContigs should be in (first ten)";
    unsigned int i = 0;
    for (auto &c: Globals::sequences) {
      std::cerr << " '" << c.first << "'";
      ++i;
      if (i >= 10) break;
    }
    std::cerr << "...\nExiting now.\n";
    exit(EXIT_FAILURE);
  }
  if (! output.empty()) {
    output.append(gap);
  }
  if (refInterval.end > contig.size()) {
    std::cerr << "Error!  Contig '" << refInterval.ref << "' has size " << contig.size() << ", but sub-sequence " << refInterval.start << "-" << refInterval.end << " is requested.\n";
    exit(EXIT_FAILURE);
  }
  // C++ is 0-based
  contig = contig.substr(refInterval.start - 1, refInterval.end - refInterval.start + 1);
  if (! refInterval.strand) {
    complement(contig);
  }
  output.append(contig);
}

void scaffolds_to_fasta (RefIntervalsSet &refIntervalsSet) {
  std::cerr << "Writing scaffolds...\n";
  std::string gap(Globals::filler_size, 'N');
  std::string sequence;
  std::ofstream output_file (Globals::output_file_name, std::ofstream::out);
  if (! output_file.is_open()){
      std::cerr << "Error!  Cannot open file '" << Globals::output_file_name << "'.\n";
      exit(EXIT_FAILURE);
  }
  for (unsigned int n_scaffolds = 0; n_scaffolds < refIntervalsSet.size(); ++n_scaffolds) {
    RefIntervals &refIntervals = refIntervalsSet[n_scaffolds];
    for (RefInterval &refInterval: refIntervals) {
      std::string current_sequence;
      extract_subsequence(refInterval, gap, current_sequence);
      sequence.append(current_sequence);
    }
    write_fasta_sequence(n_scaffolds + 1, sequence, output_file);
    sequence.clear();
    if (n_scaffolds % 100 == 0) {
      std::cerr << TAB << n_scaffolds << "/" << refIntervalsSet.size() << " scaffolds\r" << std::flush;
    }
  }
  output_file.close();
  std::cerr << TAB << refIntervalsSet.size() << "/" << refIntervalsSet.size() << " scaffolds\n";
}
