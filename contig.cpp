#include <unordered_set>
#include <sstream>
#include <fstream>

#include "constants.h"
#include "globals.h"
#include "contig.h"


// Compute the Jaccard index between two sets of barcodes
void intersectMoleculesSize(std::vector < unsigned long int > &b1, std::vector < unsigned long int > &b2, unsigned int &n_inter, double &jaccard){
  std::vector < unsigned long int > intersection_barcodes;
  std::vector < unsigned long int > union_barcodes;
  set_intersection(b1.begin(), b1.end(), b2.begin(), b2.end(), std::back_inserter(intersection_barcodes));
  set_union(b1.begin(), b1.end(), b2.begin(), b2.end(), std::back_inserter(union_barcodes));
  unsigned int n_union = union_barcodes.size();
  n_inter              = intersection_barcodes.size();
  jaccard              = 0.0;
  if (n_union == 0) return;
  jaccard = static_cast < double > (n_inter) / n_union;
}

// Find the number of common barcodes between contig ends
void add_molecules_to_contigs_extremites (Molecules &molecules, Contigs &contigs) {
  std::cerr << TAB << "Anchoring molecules...\n";
  unsigned long int barcode_id;
  std::unordered_map < std::string, unsigned long int > barcode_to_id;
  unsigned int n_barcodes_begin  = 0;
  unsigned int n_barcodes_end    = 0;
  unsigned int n_barcodes_other  = 0;
  unsigned int n_barcodes_unused = 0;
  std::ofstream mapping_file;
  if (! Globals::mapping_file_name.empty()) {
    mapping_file.open(Globals::mapping_file_name);
  }
  long unsigned int i;
  for (i = 0; i < molecules.size(); ++i) {
    Molecule &molecule = molecules[i];
    auto pos = barcode_to_id.find(molecule.barcode);
    if (pos == barcode_to_id.end()) {
      barcode_id = barcode_to_id.size();
      barcode_to_id[molecule.barcode] = barcode_to_id.size();
    }
    else {
      barcode_id = pos->second;
    }
    for (ContigPart &contigPart: contigs[molecule.chrid].contigParts) {
      if (contigPart.get_overlap(molecule) >= Globals::min_overlap) {
        // This is the condition to set a barcode to the beginning of a contig part
        if ((molecule.start >= contigPart.start) &&
            (molecule.start <= contigPart.start + std::min(Globals::window, contigPart.getSize() / 2)) &&
            (molecule.end <= contigPart.start + contigPart.getSize() * (1 - Globals::beginning_ratio))) {
          contigPart.add_beg_molecule(barcode_id);
          ++n_barcodes_begin;
          if (! Globals::mapping_file_name.empty()) {
            mapping_file << "Barcode " << molecule.barcode << TAB << Globals::chrs[molecule.chrid] << TAB << molecule.start << TAB << molecule.end << TAB << molecule.n_reads << " reads" << TAB << Globals::chrs[molecule.chrid] << TAB << contigPart.start << TAB << contigPart.end << TAB << "begin\n";
          }
        }
        // This is the condition to set a barcode to the end of a contig part
        else if ((molecule.end <= contigPart.end) &&
            (molecule.start >= contigPart.start + contigPart.getSize() * Globals::beginning_ratio) && 
            (molecule.end >= contigPart.end - std::min(Globals::window, contigPart.getSize() / 2))) {
          contigPart.add_end_molecule(barcode_id);
          ++n_barcodes_end;
          if (! Globals::mapping_file_name.empty()) {
            mapping_file << "Barcode " << molecule.barcode << TAB << Globals::chrs[molecule.chrid] << TAB << molecule.start << TAB << molecule.end << TAB << molecule.n_reads << " reads" << TAB << Globals::chrs[molecule.chrid] << TAB << contigPart.start << TAB << contigPart.end << TAB << "end\n";
          }
        }
        else {
          contigPart.add_other_molecule(barcode_id);
          ++n_barcodes_other;
          if (! Globals::mapping_file_name.empty()) {
            mapping_file << "Barcode " << molecule.barcode << TAB << Globals::chrs[molecule.chrid] << TAB << molecule.start << TAB << molecule.end << TAB << molecule.n_reads << " reads" << TAB << Globals::chrs[molecule.chrid] << TAB << contigPart.start << TAB << contigPart.end << TAB << "other\n";
          }
        }
      }
      else {
        ++n_barcodes_unused;
      }
    }
    if (i % 1000000 == 0) std::cerr << TAB << TAB << i << "/" << molecules.size() << " molecules\r" << std::flush;
  }
  std::cerr << TAB << TAB << molecules.size() << "/" << molecules.size() <<  " molecules\n";
  std::cerr << TAB << TAB << TAB << n_barcodes_begin << " anchored on the left part\n";
  std::cerr << TAB << TAB << TAB << n_barcodes_end << " anchored on the right part\n";
  std::cerr << TAB << TAB << TAB << (n_barcodes_begin + n_barcodes_end + n_barcodes_other) << " anchored in general\n";
  std::cerr << TAB << TAB << TAB << n_barcodes_unused << " unused\n";
  std::cerr << TAB << TAB << TAB << barcode_to_id.size() << " different barcodes seen\n";
  std::cerr << TAB << "Sorting barcodes...\n";
  for (size_t i = 0; i < contigs.size(); ++i) {
    Contig &contig = contigs[i];
    contig.sort_barcodes();
    if (i % 100 == 0) std::cerr << TAB << TAB << "Contig " << i << "/" << contigs.size() << "\r" << std::flush;
  }
  std::cerr << TAB << TAB << "Contig " << contigs.size() << "/" << contigs.size() << "\n";
}


Contig &get_contig (Contigs &contigs, NodeId &nodeId) {
  return contigs[nodeId.contigId];
}

ContigPart &get_contig_part (Contigs &contigs, NodeId &nodeId) {
  return contigs[nodeId.contigId].contigParts[nodeId.contigPartId];
}
