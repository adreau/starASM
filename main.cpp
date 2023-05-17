#include <iostream>

#include "globals.h"
#include "parse_parameters.h"
#include "read.h"
#include "sam_parser.h"
#include "parsers.h"
#include "fasta.h"
#include "create_molecules.h"
#include "splitASM.h"
#include "joinASM.h"

int main (int argc, char* argv[]) {
  parse_parameters(argc, argv);

  Molecules molecules;
  Sequences sequences;

  read_fasta();

  // start with SAM file
  if (Globals::output_molecules_file_name.empty()) {
    make_molecules(molecules);
    split(molecules);
    join(molecules);
  }
  // start with molecule file
  else if (Globals::output_split_file_name.empty()) {
    parse_molecule_file(molecules);
    split(molecules);
    join(molecules);
  }
  // start with split file
  else {
    RefIntervalsSet refIntervalsSet;
    parse_split_file(refIntervalsSet);
    join(molecules);
  }
  return EXIT_SUCCESS;
}
