#include "globals.h"
#include "read.h"
#include "sam_parser.h"
#include "mol_stat.h"
#include "outliers_det.h"
#include "create_molecules.h"
#include "parsers.h"
#include "molsplit.h"

void initialize_stats (Molecule_stats &molecule_stats) {
  molecule_stats.resize(Globals::chrs.size());
  for (size_t chrid = 0; chrid < Globals::chrs.size(); ++chrid) {
    long int size = (Globals::chr_sizes[chrid] - 1) / Globals::window + 1;
    molecule_stats[chrid].resize(size);
  }
}

void split () {
  Molecules      molecules;
  Molecule_stats molecule_stats;

  if (! Globals::output_molecules_file_name.empty()) {
    parse_contig_file();
    parse_molecule_file(molecules);
  }
  else {
    Barcodes barcodes;
    read_sam(barcodes);
    make_molecules(barcodes, molecules);
  }
  initialize_stats(molecule_stats);
  compute_stats(molecule_stats, molecules);
  detect_outliers(molecule_stats);
}
