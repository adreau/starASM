#include "globals.h"
#include "contig.h"
#include "scaffold.h"
#include "scaffolds_to_fasta.h"
#include "joinASM.h"

void join (Molecules &molecules, Contigs &contigs) {
  std::cerr << "Joining contigs...\n";
  add_molecules_to_contigs_extremites(molecules, contigs);
  Graph graph;
  create_nodes(contigs, graph);
  create_arcs(contigs, graph);
  if (! Globals::graph_file_name.empty()) {
    write_graph(contigs, graph);
  }
  remove_bifurcations(graph);
  Scaffolds scaffolds;
  find_scaffolds(graph, scaffolds);
  RefIntervalsSet refIntervalsSet;
  scaffolds_to_intervals(scaffolds, contigs, refIntervalsSet);
  merge_close_contigs(refIntervalsSet);
  if (! Globals::scaffold_file_name.empty()) {
    print_scaffold(refIntervalsSet);
  }
  scaffolds_to_fasta(refIntervalsSet);
}
