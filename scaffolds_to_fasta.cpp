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

void extract_subsequence (std::map < std::string, std::string > &contigs, RefInterval &refInterval, std::string &gap, std::string &output) {
  std::string contig;
  output.clear();
  try {
    contig = contigs.at(refInterval.ref);
  }
  catch (const std::out_of_range& oor) {
    std::cerr << "Cannot find contig '" << refInterval.ref << "'.\nContigs should be in (first ten)";
    unsigned int i = 0;
    for (auto &c: contigs) {
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
  std::ofstream scaff_fasta(Globals::output_file_name, std::ofstream::out);	
  if (! scaff_fasta.is_open()) {
    std::cerr << "Error!  Cannot open file '" << Globals::output_file_name << "'" << std::endl;
    exit(EXIT_FAILURE);
  }

  std::map < std::string, std::string > contigs;
  read_fasta(contigs, Globals::contigs_file_name);

  std::string gap(Globals::filler_size, 'N');
  std::string sequence;
  unsigned int n_scaffolds;

  for (n_scaffolds = 0; n_scaffolds < refIntervalsSet.size(); ++n_scaffolds) {
    RefIntervals &refIntervals = refIntervalsSet[n_scaffolds];
    for (RefInterval &refInterval: refIntervals) {
      std::string current_sequence;
      extract_subsequence(contigs, refInterval, gap, current_sequence);
      sequence.append(current_sequence);
    }
    write_fasta_sequence(n_scaffolds + 1, sequence, scaff_fasta);
    sequence.clear();

    if (n_scaffolds % 100 == 0) {
      std::cout << n_scaffolds << " lines read.\r" << std::flush;
    }
  }
  std::cout << n_scaffolds << " lines read, done.\n";
}
