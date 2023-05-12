#include <cassert>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

#include "globals.h"
#include "interval.h"
#include "parsers.h"

void parse_contig_parts (RefIntervalsSet &refIntervalsSet) {
  std::ifstream joins_file (Globals::joins_file_name);
  if (! joins_file.is_open()) {
    std::cerr << "Error!  Cannot open file '" << Globals::joins_file_name << "'" << std::endl;
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
