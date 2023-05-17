#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cassert>

#include "globals.h"
#include "parse_parameters.h"
#include "contig.h"
#include "graph.h"
#include "scaffold.h"
#include "scaffolds_to_fasta.h"
#include "joinASM.h"

void join (Molecules &molecules, Contigs &contigs) {
  add_molecules_to_contigs_extremites(molecules, contigs);
  Graph graph;
  create_nodes(contigs, graph);
  create_arcs(contigs, graph);
  write_graph(contigs, graph);
  remove_bifurcations(graph);
  Scaffolds scaffolds;
  find_scaffolds(graph, scaffolds);
  RefIntervalsSet refIntervalsSet;
  scaffolds_to_intervals(scaffolds, contigs, refIntervalsSet);
  merge_close_contigs(refIntervalsSet);
  print_scaffold(refIntervalsSet);
  scaffolds_to_fasta(refIntervalsSet);
}
