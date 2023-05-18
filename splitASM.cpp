#include "globals.h"
#include "read.h"
#include "sam_parser.h"
#include "mol_stat.h"
#include "outliers_det.h"
#include "create_molecules.h"
#include "parsers.h"
#include "splitASM.h"

void initialize_stats (Molecule_stats &molecule_stats) {
  molecule_stats.resize(Globals::chrs.size());
  for (size_t chrid = 0; chrid < Globals::chrs.size(); ++chrid) {
    long int size = (Globals::chr_sizes[chrid] - 1) / Globals::window + 1;
    molecule_stats[chrid].resize(size);
  }
}

void split (Molecules &molecules, Contigs &contigs) {
  std::cerr << "Splitting contigs...\n";
  Molecule_stats molecule_stats;
  initialize_stats(molecule_stats);
  std::cerr << TAB << "Computing stats...\n";
  compute_stats(molecule_stats, molecules);
  detect_outliers(molecule_stats, contigs);
}
