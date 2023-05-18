#include <fstream>

#include "interval.h"
#include "graph.h"


// Merge two consecutive contigs if they belong to the same scaffold, have been previously split, and are not too distant from each other
void merge_close_contigs (RefIntervalsSet &refIntervalsSet) {
  unsigned int n_merges = 0;
  for (RefIntervals &refIntervals: refIntervalsSet) {
    size_t i = 0;
    for (size_t j = 1; j < refIntervals.size(); ++j) {
      if (refIntervals[i].can_merge(refIntervals[j])) {
        refIntervals[i].merge(refIntervals[j]);
        ++n_merges;
      }
      else {
        ++i;
        if (i != j) {
          refIntervals[i] = refIntervals[j];
        }
      }
    }
    if (i < refIntervals.size()) {
      refIntervals.erase(refIntervals.begin() + i + 1, refIntervals.end());
    }
  }
  std::cerr << n_merges << " contig part merges.\n";
}


void print_scaffold (RefIntervalsSet &refIntervalsSet) {
  std::ofstream scaffold_file (Globals::scaffold_file_name, std::ofstream::out);
  if (! scaffold_file.is_open()){
      std::cerr << "Error!  Cannot open file '" << Globals::scaffold_file_name << "'.\n";
      exit(EXIT_FAILURE);
  }
  for (auto &refIntervals: refIntervalsSet) {
    for (auto &refInterval: refIntervals) {
      scaffold_file << refInterval << ';';
    }
    scaffold_file << '\n';
  }
  scaffold_file.close();
}
