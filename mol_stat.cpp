#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include "mol_stat.h"

void compute_stats (Molecule_stats &molecule_stats, Molecules &molecules) {
  unsigned long cpt = 0;
  for (auto &molecule: molecules) {
    // Molecules are 0-based
    unsigned int window_start = molecule.start / Globals::window;
    unsigned int window_end   = molecule.end   / Globals::window;
    unsigned int molecule_size = molecule.end - molecule.start + 1;
    assert(molecule.chrid < molecule_stats.size());
    assert(window_end < molecule_stats[molecule.chrid].size());
    for (unsigned int windowid = window_start; windowid <= window_end; ++windowid) {
      unsigned int beg_window = std::max < unsigned int > (windowid * Globals::window, molecule.start);
      unsigned int end_window = std::min < unsigned int > (beg_window + Globals::window - 1, molecule.end);
      unsigned int size = end_window - beg_window + 1;
      molecule_stats[molecule.chrid][windowid].coverage     += size;
      molecule_stats[molecule.chrid][windowid].length       += molecule_size * size;
      molecule_stats[molecule.chrid][windowid].read_density += static_cast < double > (molecule.n_reads) / molecule_size * size;
    }
    ++molecule_stats[molecule.chrid][window_start].start;
    ++molecule_stats[molecule.chrid][window_end].end;
    ++cpt;
    if (cpt % 10000000 == 0) std::cerr << TAB << TAB << cpt << "/" << molecules.size() << " molecules\r" << std::flush;
  }
  std::cerr << TAB << TAB << molecules.size() << "/" << molecules.size() << " molecules\n";
  std::ofstream counts_file;
  if (! Globals::counts_file_name.empty()) counts_file.open(Globals::counts_file_name, std::ofstream::out);
  for (size_t chrid = 0; chrid < Globals::chrs.size(); ++chrid) {
    std::string &chr = Globals::chrs[chrid];
    unsigned long size = Globals::chr_sizes[chrid] / Globals::window;
    for (size_t windowid = 0; windowid <= size; ++windowid) {
      unsigned int beg_window = windowid * Globals::window;
      unsigned int end_window = std::min < int > (beg_window + Globals::window - 1, Globals::chr_sizes[chrid]);
      double size_dbl = static_cast < double > (end_window - beg_window + 1);
      // Careful in the order of normalization
      molecule_stats[chrid][windowid].coverage   /= size_dbl;
      if (molecule_stats[chrid][windowid].coverage == 0) {
        molecule_stats[chrid][windowid].length       = 0;
        molecule_stats[chrid][windowid].read_density = 0;
      }
      else {
        molecule_stats[chrid][windowid].length       /= molecule_stats[chrid][windowid].coverage * Globals::window;
        molecule_stats[chrid][windowid].read_density /= molecule_stats[chrid][windowid].coverage * Globals::window;
      }

      // Switch back to 1-based positions
      if (! Globals::counts_file_name.empty()) {
        counts_file << chr                             << TAB <<
          windowid * Globals::window + 1               << TAB <<
          (windowid + 1) * Globals::window             << TAB <<
          molecule_stats[chrid][windowid].coverage     << TAB <<
          molecule_stats[chrid][windowid].length       << TAB <<
          molecule_stats[chrid][windowid].read_density << TAB <<
          molecule_stats[chrid][windowid].start        << TAB <<
          molecule_stats[chrid][windowid].end          << '\n';
      }
    }
  }
  if (! Globals::counts_file_name.empty()) counts_file.close();
}
