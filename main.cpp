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

  Sequences sequences;
  read_fasta();

  // start with BAM file
  if (Globals::input_molecules_file_name.empty()) {
    Molecules molecules;
    Contigs contigs;
    make_molecules(molecules);
    split(molecules, contigs);
    join(molecules, contigs);
  }
  // start with molecule file
  else if (Globals::input_split_file_name.empty()) {
    Molecules molecules;
    Contigs contigs;
    parse_molecule_file(molecules);
    split(molecules, contigs);
    join(molecules, contigs);
  }
  // start with split file
  else {
    Molecules molecules;
    Contigs contigs;
    parse_molecule_file(molecules);
    parse_split_file(contigs);
    join(molecules, contigs);
  }
  std::cerr << "Done.\n";
  return EXIT_SUCCESS;
}
