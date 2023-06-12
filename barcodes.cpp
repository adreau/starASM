#include <algorithm>
#include <iostream>
#include "constants.h"
#include "globals.h"
#include "barcodes.h"

unsigned int Barcodes::count_n_reads (unsigned long id) {
  std::vector < unsigned long > names;
  for (unsigned i = offsets[id]; i < offsets[id] + counts[id]; ++i) {
    names.push_back(reads[i].name);
  }
  return count_different(names);
}

// Set to zero reads molecules with few reads
void Barcodes::trim () {
  std::cerr << TAB << "Removing outlier barcodes...\n";
  unsigned int cpt = 0;
  for (unsigned i = 0; i < ids.size(); ++i) {
    if (count_n_reads(i) < Globals::min_n_reads_barcode) {
      counts[i] = 0;
      ++cpt;
    }
  }
  std::cerr << TAB << TAB << cpt << "/" << ids.size() << " barcodes removed.\n";
}

void Barcodes::sort () {
  std::cerr << TAB << "Sorting barcodes...\n";
  for (unsigned i = 0; i < ids.size(); ++i) {
    if (counts[i] > 1) {
      std::sort(reads.begin() + offsets[i], reads.begin() + offsets[i] + counts[i]);
    }
  }
}

void Barcodes::set_structure () {
  unsigned long size   = count_map.size();
  unsigned int cpt     = 0;
  unsigned long offset = 0;
  offsets.reserve(size);
  counts.reserve(size);
  for (auto &it: count_map) {
    offsets.push_back(offset);
    counts.push_back(it.second);
    offset       += it.second;
    ids[it.first] = cpt;
    ++cpt;
  }
  count_map.clear();
  current_offsets = offsets;
  reads = std::vector < Read > (offset);
}
