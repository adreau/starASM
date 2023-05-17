#include "scaffold.h"

// Find an element which has at most one neighbor
size_t find_first_scaffold_end (Graph &graph, size_t nodeIdStart, std::vector < bool > &seen_nodes) {

  for (size_t nodeId = nodeIdStart; nodeId < graph.nodes.size(); ++nodeId) {
    if ((! seen_nodes[nodeId]) && (graph.nodes[nodeId].get_n_edges() < 2)) {
      return nodeId;
    }
  }
  return graph.nodes.size();
}


// Find the list of scaffolds in the graph
void find_scaffolds (Graph &graph, Scaffolds &scaffolds) {

  size_t n_nodes = graph.nodes.size();
  std::vector < bool > seen_nodes (n_nodes, false);
  for (size_t nodeIdStart = find_first_scaffold_end(graph, 0, seen_nodes); nodeIdStart < n_nodes; nodeIdStart = find_first_scaffold_end(graph, ++nodeIdStart, seen_nodes)) {
    size_t prev_node_id = nodeIdStart;
    Node &node = graph.nodes[nodeIdStart];
    seen_nodes[nodeIdStart] = true;
    scaffolds.emplace_back();
    // First of a chain
    if (node.get_n_edges() > 0) {
      // Find whether we should go right, or reverse the first contig
      Edge &edge = node.get_first_edge();
      bool is_forward = ((edge.link_type == Link_types::EB) || (edge.link_type == Link_types::EE));
      scaffolds.back().emplace_back(graph.nodes[nodeIdStart], is_forward);
      while (true) {
        size_t nodeId = edge.nodeId;
        Node &node = graph.nodes[nodeId];
        // If the edge is forward, the next edge is forward iff the link starts with E.
        // If the edge is reverse, the next edge is reverse iff the link starts with B.
        is_forward = ((edge.link_type == Link_types::EB) || (edge.link_type == Link_types::EE));
        scaffolds.back().emplace_back(node, is_forward);
        seen_nodes[nodeId] = true;
        // We reached the end of the chain
        if (node.get_n_edges() == 1) {
          break;
        }
        edge = node.get_other_edge(prev_node_id, reverse_link_type[edge.link_type]);
        // This only should happen if the scaffold is circular
        if (seen_nodes[edge.nodeId]) {
          break;
        }
        prev_node_id = nodeId;
      }
    }
    // Lone node
    else {
      scaffolds.back().emplace_back(graph.nodes[nodeIdStart], true);
    }
  }
}


void scaffolds_to_intervals (Scaffolds &scaffolds, Contigs &contigs, RefIntervalsSet &refIntervalsSet) {
  for (Scaffold &scaffold: scaffolds) {
    RefIntervals refIntervals;
    for (ScaffoldPart &scaffold_part: scaffold) {
      Contig     &contig      = get_contig(contigs, scaffold_part.nodeId);
      ContigPart &contig_part = get_contig_part(contigs, scaffold_part.nodeId);
      if (contig_part.is_set()) {
        refIntervals.emplace_back(contig.name, scaffold_part.is_forward, contig_part.start, contig_part.end);
      }
    }
    refIntervalsSet.push_back(refIntervals);
  }
}

