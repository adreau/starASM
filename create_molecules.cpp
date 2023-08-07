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
  for (auto &it: barcodes.ids) {
    unsigned long i = it.second;
    if (barcodes.counts[i] > 0) {
      unsigned long start = barcodes.offsets[i];
      unsigned long end   = barcodes.offsets[i] + barcodes.counts[i];
	  assert(end <= barcodes.reads.size());
      do {
        unsigned int chrid   = barcodes.reads[start].chrid;
        unsigned int chr_end = find_chr_end(barcodes, chrid, start, end);
		assert(chr_end <= barcodes.reads.size());
        _join_to_molecules(molecules, barcodes, it.first, chrid, start, chr_end);
        start = chr_end;
	    assert(start <= barcodes.reads.size());
      }
      while (start < end);
    }
  }
}

void remove_low_barcode_counts (Molecules &molecules) {
  std::cerr << TAB << "Removing outlier molecules.\n";
  // Compute the number of times each barcode is seen
  std::unordered_map < std::string, unsigned int > counts_map;
  for (Molecule &molecule: molecules) {
    ++counts_map[molecule.barcode];
  }
  unsigned int max_count = 0;
  for (auto &it: counts_map) {
    max_count = std::max(max_count, it.second);
  }
  // Compute the distribution of the occurrences
  std::vector < unsigned int > n_occurrences (max_count + 1);
  for (auto &it: counts_map) {
    ++n_occurrences[it.second];
  }
  unsigned int min_threshold; // Values before this threshold should be removed
  unsigned int max_threshold; // Actually, the maximum of the "reliable" distribution, times 3
  // Find the first decrease of the distribution
  unsigned int previous_value = 0;
  unsigned long id = 0;
  for (; id < n_occurrences.size(); ++id) {
std::cerr << "step 1: " << id << " => " << n_occurrences[id] << "/" << n_occurrences.size() << " vs " << previous_value << "\n";
    if (n_occurrences[id] < previous_value) {
	  break;
    }
    previous_value = n_occurrences[id];
  }
  if (id == n_occurrences.size()) {
    std::cerr << TAB << TAB << "Cannot model the distribution (step 1).  Skipping this step.\n";
	return;
  }
  // Find the next increase of the distribution
  for (; id < n_occurrences.size(); ++id) {
std::cerr << "step 2: " << id << " => " << n_occurrences[id] << "/" << n_occurrences.size() << " vs " << previous_value << "\n";
    if (n_occurrences[id] > previous_value) {
	  break;
    }
    previous_value = n_occurrences[id];
  }
  if (id == n_occurrences.size()) {
    std::cerr << TAB << TAB << "Cannot model the distribution (step 2).  Skipping this step.\n";
	return;
  }
  min_threshold = id;
  // Find the next decrease of the distribution
  for (; id < n_occurrences.size(); ++id) {
std::cerr << "step 3: " << id << " => " << n_occurrences[id] << "/" << n_occurrences.size() << " vs " << previous_value << "\n";
    if (n_occurrences[id] < previous_value) {
	  break;
    }
    previous_value = n_occurrences[id];
  }
  if (id == n_occurrences.size()) {
    std::cerr << TAB << TAB << "Cannot model the distribution (step 3).  Skipping this step.\n";
	return;
  }
  max_threshold = id;
  max_threshold *= 3;
  if (min_threshold > 20) {
    std::cerr << TAB << TAB << "No barcode threshold applied.\n";
	return;
  }
  // In-place trimming
  unsigned int n_barcodes_removed = 0;
  for (auto &it: counts_map) {
    if ((it.second < min_threshold) || (it.second > max_threshold)) {
      ++n_barcodes_removed;
	}
  }
  unsigned int n_molecules_removed = 0;
  unsigned long write_molecule_id  = 0;
  for (unsigned long read_molecule_id = 0; read_molecule_id < molecules.size(); ++read_molecule_id) {
	unsigned long count = counts_map[molecules[read_molecule_id].barcode];
    if ((count < min_threshold) || (count > max_threshold)) {
      ++n_molecules_removed;
	}
	else {
      if (write_molecule_id != read_molecule_id) {
        molecules[write_molecule_id] = molecules[read_molecule_id];
	  }
	  ++write_molecule_id;
	}
  }
  std::cerr << TAB << TAB << "Barcode thresholds: " << min_threshold << "--" << max_threshold << ", removed " << n_barcodes_removed << "/" << counts_map.size() << " barcodes and " << n_molecules_removed << "/" << molecules.size() << " molecules.\n";
  molecules.resize(write_molecule_id + 1);
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
  barcodes.remove_low_read_counts();
  barcodes.sort();
  join_to_molecules(barcodes, molecules);
  remove_low_barcode_counts(molecules);
  sort_molecules(molecules);
  if (! Globals::output_molecules_file_name.empty()) {
    print_molecules(molecules);
  }
}
