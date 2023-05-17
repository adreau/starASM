#ifndef SCAFFOLD_H
#define SCAFFOLD_H

#include "interval.h"
#include "graph.h"

// A scaffold part is an id to a contig part, and a strand (reversed, or not)
struct ScaffoldPart {
  NodeId nodeId;
  bool   is_forward;

  ScaffoldPart (NodeId &n, bool b): nodeId(n), is_forward(b) {}

  friend std::ostream &operator<< (std::ostream &out, ScaffoldPart &sp);
};

inline std::ostream &operator<< (std::ostream &out, ScaffoldPart &sp) {
  out << (sp.is_forward? '+': '-') << sp.nodeId;
  return out;
}

using Scaffold = std::vector < ScaffoldPart >;

using Scaffolds = std::vector < Scaffold >;

void find_scaffolds (Graph &graph, Scaffolds &scaffolds);
void scaffolds_to_intervals (Scaffolds &scaffolds, Contigs &contigs, RefIntervalsSet &refIntervalsSet);

#endif
