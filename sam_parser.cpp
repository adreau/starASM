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
bool parse_main_line (std::string &line, std::string &name, unsigned int &chrid, unsigned long &start, unsigned long &end, unsigned int &mapq, std::string &barcode) {
  std::stringstream split_line(line);
  std::string value;
  const unsigned int no_value = -1;
  unsigned int flag = no_value;
  for (unsigned int i = 0; split_line >> value; ++i) {
    if (i == 0) {
      name = value;
    }
    else if (i == 1) {
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

void add_barcode (Barcodes &barcodes, std::string &name, unsigned int chrid, unsigned long start, unsigned long end, unsigned int mapq, std::string &barcode) {
  if (barcodes.find(barcode) == barcodes.end()) {
    barcodes[barcode] = std::vector < std::vector < Read > > (Globals::chrs.size());
  }
  assert(chrid < Globals::chrs.size());
  barcodes[barcode][chrid].emplace_back(name, start, end, mapq);
}

// Returns true iff read passed the filters
bool read_main_line (Barcodes &barcodes, std::string &line) {
  std::string name;
  unsigned int chrid = -1;
  unsigned long start = -1, end = -1;
  unsigned int mapq = -1;
  std::string barcode;
  if (line[0] == '@') {
    return false;
  }
  if (parse_main_line(line, name, chrid, start, end, mapq, barcode)) {
    add_barcode(barcodes, name, chrid, start, end, mapq, barcode);
    return true;
  }
  return false;
}

void read_sam (Barcodes &barcodes) {
  if (Globals::input_file_name.empty()) {
    std::cerr << "Error!  Input SAM file is missing.\n";
    exit(EXIT_FAILURE);
  }
  std::ifstream input_file(Globals::input_file_name);
  if (! input_file.is_open()) {
    std::cerr << "Error!  Cannot open file SAM file '" << Globals::input_file_name << "'.\n";
    exit(EXIT_FAILURE);
  }
  unsigned long int n_reads_kept = 0;
  unsigned long int cpt = 0;
  std::cerr << "Reading SAM file...\n";
  std::string line;
  for (; std::getline(input_file, line); ++cpt) {
    if (read_main_line(barcodes, line)) {
      ++n_reads_kept;
    }
    ++cpt;
    if (cpt % 10000000 == 0) {
      std::cerr << TAB << cpt << " lines read, " << n_reads_kept << " reads kept, using " << barcodes.size() << " barcodes.\r" << std::flush;
    }
  }
  std::cerr << TAB << cpt << " lines read, " << n_reads_kept << " reads kept, using " << barcodes.size() << " barcodes.\n";
}
