#ifndef BARCODES_H
#define BARCODES_H

#include <string>
#include <unordered_map>
#include "read.h"


struct Barcodes {
  std::unordered_map < std::string, unsigned int > count_map;
  std::unordered_map < std::string, unsigned int > ids;
  std::vector < unsigned int >                     counts;
  std::vector < unsigned int >                     offsets;
  std::vector < unsigned int >                     current_offsets;
  std::vector < Read >                             reads;

  void          set_structure();
  unsigned int  count_n_reads (unsigned long id);
  void          remove_low_read_counts();
  void          sort();
  void          remove_low_barcode_counts();
};

#endif
