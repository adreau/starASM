#include <cassert>
#include <algorithm>    // std::sort
#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>     // std::cin, std::cout, std::cerr
#include <sstream>      // std::stringstream
#include <fstream>

#include "constants.h"
#include "globals.h"
#include "parse_parameters.h"
#include "sam_parser.h"
#include "create_molecules.h"



/*
unsigned int count_n_names (std::vector < unsigned long > &names) {
  std::sort(names.begin(), names.end());
  auto last = std::unique(names.begin(), names.end());
  return std::distance(names.begin(), last);
}

unsigned int count_n_reads (std::vector < std::vector < Read > > &reads) {
  std::vector < unsigned long > names;
  for (auto &r1: reads) {
    for (auto &r2: r1) {
      names.push_back(r2.name);
    }
  }
  return count_n_names(names);
}

void trim_barcodes (Barcodes &barcodes) {
  std::cerr << TAB << "Removing outlier barcodes...\n";
  auto it = barcodes.begin();
  while (it != barcodes.end()) {
    if (count_n_reads(it->second) < Globals::min_n_reads_barcode) it = barcodes.erase(it);
    else                                                        ++it;
  }
}

void sort_barcodes (Barcodes &barcodes) {
  std::cerr << TAB << "Sorting " << barcodes.size() << " barcodes...\n";
  for (auto &p: barcodes) {
    p.second.shrink_to_fit();
    for (auto &q: p.second) {
      q.shrink_to_fit();
      std::sort(q.begin(), q.end());
    }
  }
}
*/

// Returns the index of the last read of the current split.
//   A split is made iff the distance between consecutive reads is greater than max_read_distance
unsigned int find_first_split (Barcodes &barcodes, unsigned int start_id, unsigned int end_id) {
  for (unsigned int prev = start_id, next = start_id + 1; next < end_id; ++prev, ++next) {
    if (barcodes.reads[prev].get_distance(barcodes.reads[next]) > Globals::max_read_distance) {
      return next;
    }
  }
  return end_id;
}

// Returns true iff at least one reads has the given mapq
bool check_min_mapq (std::vector < Read > &reads, unsigned int start_id, unsigned int end_id) {
  assert(start_id <= end_id);
  assert(end_id <= reads.size());
  unsigned int mapq = 0;
  for (unsigned int id = start_id; id < end_id; ++id) {
    mapq = std::max(mapq, reads[id].mapq);
  }
  return (mapq >= Globals::min_mapq_solid);
}

// Returns true iff at least one solid read is found
bool find_solid_ends (std::vector < Read > &reads, unsigned int start_id, unsigned int end_id, unsigned int &start_solid_id, unsigned int &end_solid_id) {
  assert(start_id <= end_id);
  assert(end_id <= reads.size());
  constexpr unsigned int no_id = -1;
  start_solid_id = end_solid_id = no_id;
  for (unsigned int id = start_id; id < end_id; ++id) {
    if (reads[id].is_solid()) {
      if (start_solid_id == no_id) {
        start_solid_id = id;
      }
      end_solid_id = id + 1;
    }
  }
  return (start_solid_id != no_id);
}

void add_molecule (Molecules &molecules, Barcodes &barcodes, unsigned int start_id, unsigned int end_id, const std::string &barcode, unsigned int chrid) {
  assert(start_id < end_id);
  assert(end_id <= barcodes.reads.size());
  std::vector < unsigned long > names;
  names.reserve(end_id - start_id);
  for (unsigned int id = start_id; id < end_id; ++id) {
    names.push_back(barcodes.reads[id].name);
  }
  assert(barcodes.reads[start_id].start < barcodes.reads[end_id - 1].end);
  molecules.emplace_back(chrid, barcodes.reads[start_id].start, barcodes.reads[end_id - 1].end, barcode, count_different(names));
}

void _join_to_molecules (Molecules &molecules, Barcodes &barcodes, const std::string &barcode, unsigned int chrid, unsigned long start_id, unsigned long end_id) {
  unsigned int next_id;
  unsigned int start_solid_id, end_solid_id;
  do {
    next_id = find_first_split(barcodes, start_id, end_id);
    assert(next_id <= barcodes.reads.size());
    if (check_min_mapq(barcodes.reads, start_id, next_id) && find_solid_ends(barcodes.reads, start_id, next_id, start_solid_id, end_solid_id)) {
      add_molecule(molecules, barcodes, start_solid_id, end_solid_id, barcode, chrid);
    }
    start_id = next_id;
  }
  while (start_id < end_id);
}

// Suppose that the reads are sorted by chrid
// Return the index of the first read that do not belong to this chr
unsigned long find_chr_end (Barcodes &barcodes, unsigned int chrid, unsigned long start, unsigned long end) {
  for (unsigned long i = start; i < end; ++i) {
    if (barcodes.reads[i].chrid != chrid) {
      return i;
    }
  }
  return end;
}

void join_to_molecules (Barcodes &barcodes, Molecules &molecules) {
  unsigned long i = 0;
  for (auto &it: barcodes.ids) {
    if (barcodes.counts[i] > 0) {
      unsigned long start = barcodes.offsets[i];
      unsigned long end   = barcodes.offsets[i] + barcodes.counts[i];
      do {
        unsigned int chrid   = barcodes.reads[start].chrid;
        unsigned int chr_end = find_chr_end(barcodes, chrid, start, end);
        _join_to_molecules(molecules, barcodes, it.first, chrid, start, chr_end);
        start = chr_end;
      }
      while (start < end);
    }
    ++i;
  }
}

void sort_molecules (Molecules &molecules) {
  std::cerr << TAB << "Sorting " << molecules.size() << " molecules...\n";
  molecules.shrink_to_fit();
  std::sort(molecules.begin(), molecules.end());
}

void print_molecules (Molecules &molecules) {
  std::ofstream output_file(Globals::output_molecules_file_name);
  if (! output_file) {
    std::cerr << "Cannot write output molecule file to " << Globals::output_molecules_file_name << ".\nExting." << "\n";
    exit(EXIT_FAILURE);
  }
  for (Molecule &molecule: molecules) {
    output_file << molecule;
  }
  output_file.close();
}

void make_molecules(Molecules &molecules) {
  Barcodes barcodes;
  read_sam(barcodes);
  std::cerr << "Creating molecules from reads...\n";
  barcodes.trim();
  barcodes.sort();
  join_to_molecules(barcodes, molecules);
  sort_molecules(molecules);
  if (! Globals::output_molecules_file_name.empty()) {
    print_molecules(molecules);
  }
}
