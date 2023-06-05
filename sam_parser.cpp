#include <cassert>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include "constants.h"
#include "globals.h"
#include "read.h"
#include "molecule.h"
#include "sam_parser.h"

bool is_mapped (unsigned int flag) {
  return ((flag & 4) == 0);
}

bool is_duplicate (unsigned int flag) {
  return ((flag & 1024) != 0);
}

// Return the match length
unsigned int parse_cigar (std::string &cigar) {
  unsigned int match_len = 0;
  unsigned int number = 0;
  for (char c: cigar) {
    switch (c) {
      case '0':
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
        number = number * 10 + (c - '0');
        break;
      case 'M':
      case 'X':
      case '=':
      case 'D':
        match_len += number;
        number = 0;
        break;
      case 'I':
      case 'S':
      case 'P':
      case 'H':
        number = 0;
        break;
      default:
        assert(false);
    }
  }
  return match_len;
}

// Returns true iff:
//  - read is mapped
//  - flagged in proper pair
//  - mapq is higher than threshold
//  - barcode in BX tag is set
bool parse_main_line (std::string &line, unsigned int &chrid, unsigned long &start, unsigned long &end, unsigned int &mapq, std::string &barcode) {
  std::stringstream split_line(line);
  std::string value;
  const unsigned int no_value = -1;
  unsigned int flag = no_value;
  for (unsigned int i = 0; split_line >> value; ++i) {
    if (i == 1) {
      flag = std::stoi(value);
      if ((! is_mapped(flag)) || (is_duplicate(flag))) {
        return false;
      }
    }
    else if (i == 2) {
      if (value == "*") {
        return false;
      }
      assert(Globals::chrids.find(value) != Globals::chrids.end());
      chrid = Globals::chrids[value];
    }
    else if (i == 3) {
      start = std::stoul(value) - 1;
    }
    else if (i == 4) {
      mapq = std::stoi(value);
      if (mapq < Globals::min_mapq) {
        return false;
      }
    }
    else if (i == 5) {
      end = start + parse_cigar(value);
    }
    else if (i >= 11) {
      if (starts_with(value, BARCODE_FLAG)) {
        barcode = value.substr(BARCODE_FLAG.size());
      }
    }
  }
  if (barcode.empty()) {
    return false;
  }
  assert(chrid != no_value);
  assert(start != no_value);
  assert(end   != no_value);
  assert(mapq  != no_value);
  return true;
}

// Add the barcode count
// Returns true iff read passed the filters
bool add_barcode_count (Barcodes &barcodes, std::string &line) {
  unsigned int chrid = -1;
  unsigned long start = -1, end = -1;
  unsigned int mapq = -1;
  std::string barcode;
  if (line[0] == '@') {
    return false;
  }
  if (parse_main_line(line, chrid, start, end, mapq, barcode)) {
std::cerr << "previous count: " << barcodes.count_map[barcode] << ", ";
    ++barcodes.count_map[barcode];
std::cerr << "now: " << barcodes.count_map[barcode] << "\n";
    return true;
  }
  return false;
}

// Count the number of occurrences of each barcode
void count_barcodes (Barcodes &barcodes) {
  std::ifstream input_file(Globals::input_file_name);
  if (! input_file.is_open()) {
    std::cerr << "Error!  Cannot open file SAM file '" << Globals::input_file_name << "'.\n";
    exit(EXIT_FAILURE);
  }
  unsigned long int n_reads_kept = 0;
  unsigned long int cpt = 0;
  std::cerr << "Reading SAM file for barcode counting...\n";
  std::string line;
  for (; std::getline(input_file, line); ++cpt) {
    if (add_barcode_count(barcodes, line)) {
      ++n_reads_kept;
    }
    ++cpt;
    if (cpt % 10000000 == 0) {
      std::cerr << TAB << cpt << " lines read, " << n_reads_kept << " reads kept, using " << barcodes.count_map.size() << " barcodes.\r" << std::flush;
    }
  }
  std::cerr << TAB << cpt << " lines read, " << n_reads_kept << " reads kept, using " << barcodes.count_map.size() << " barcodes.\n";
}

void add_barcode (Barcodes &barcodes, unsigned int chrid, unsigned long start, unsigned long end, unsigned int mapq, const std::string &barcode) {
  auto it = barcodes.ids.find(barcode);
  if (it == barcodes.ids.end()) {
    return;
  }
  unsigned int id = it->second;
  barcodes.reads[barcodes.current_offsets[id]] = Read(1000, chrid, start, end, mapq);
  ++barcodes.current_offsets[id];
  assert(chrid < Globals::chrs.size());
}

void add_barcode_line (Barcodes &barcodes, std::string &line) {
  unsigned int chrid = -1;
  unsigned long start = -1, end = -1;
  unsigned int mapq = -1;
  std::string barcode;
  if (line[0] == '@') {
    return;
  }
  if (parse_main_line(line, chrid, start, end, mapq, barcode)) {
    add_barcode(barcodes, chrid, start, end, mapq, barcode);
  }
}

void add_barcodes (Barcodes &barcodes) {
  std::ifstream input_file(Globals::input_file_name);
  if (! input_file.is_open()) {
    std::cerr << "Error!  Cannot open file SAM file '" << Globals::input_file_name << "'.\n";
    exit(EXIT_FAILURE);
  }
  unsigned long int cpt = 0;
  std::cerr << "Reading SAM file for barcode storing...\n";
  std::string line;
  for (; std::getline(input_file, line); ++cpt) {
    add_barcode_line(barcodes, line);
    ++cpt;
    if (cpt % 10000000 == 0) {
      std::cerr << TAB << cpt << " lines read.\r" << std::flush;
    }
  }
  barcodes.current_offsets.clear();
  std::cerr << TAB << cpt << " lines read.\n";
}

void read_sam (Barcodes &barcodes) {
  if (Globals::input_file_name.empty()) {
    std::cerr << "Error!  Input SAM file is missing.\n";
    exit(EXIT_FAILURE);
  }
  count_barcodes(barcodes);
  barcodes.set_structure();
  add_barcodes(barcodes);
}
