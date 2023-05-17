#include <unordered_set>
#include <sstream>
#include <fstream>

#include "globals.h"
#include "contig.h"


// Count the number of common barcodes between to sets
unsigned int intersectMoleculesSize(std::vector < unsigned long int > &b1, std::vector < unsigned long int > &b2){

  std::vector < unsigned long int > common_barcodes;

  set_intersection(b1.begin(), b1.end(),
          b2.begin(), b2.end(),
          std::back_inserter(common_barcodes));
  size_t s1 = b1.size();
  size_t s2 = b2.size();
  size_t sc = common_barcodes.size();
  if (sc == 0) return 0;

  switch (Globals::condition) {
      case 1:
          if ((sc >= s1 * 0.8) && (sc >= s2 * 0.8)) return sc;
          return 0;
      case 2:
          if ((sc >= s1 * 0.8) || (sc >= s2 * 0.8)) return sc;
          return 0;
      case 3:
          if ((sc >= s1 * 0.6) && (sc >= s2 * 0.6)) return sc;
          return 0;
      case 4:
          if ((sc >= s1 * 0.6) || (sc >= s2 * 0.6)) return sc;
          return 0;
      case 5:
          if ((sc >= s1 * 0.4) && (sc >= s2 * 0.4)) return sc;
          return 0;
      case 6:
          if ((sc >= s1 * 0.4) || (sc >= s2 * 0.4)) return sc;
          return 0;
      case 7:
          if ((sc >= s1 * 0.2) && (sc >= s2 * 0.2)) return sc;
          return 0;
      case 8:
          if ((sc >= s1 * 0.2) || (sc >= s2 * 0.2)) return sc;
          return 0;
  }
  std::cerr << "Error: arc condition should be between 1 and 8.\n";
  exit(EXIT_FAILURE);
  return 0;
}

// Find the number of common barcodes between contig ends
void add_molecules_to_contigs_extremites (Contigs &contigs) {

  std::cerr << "Reading molecule file...\n";
  std::ifstream molecule_file (Globals::molecule_file_name);
  std::string molecule_line, barcode, ctg, prevCtg;
  unsigned long int beg_pos, end_pos, nReads, barcode_id;
  std::unordered_map < std::string, unsigned long int > barcode_to_id;
  unsigned int n_barcodes_begin  = 0;
  unsigned int n_barcodes_end    = 0;
  unsigned int n_barcodes_other  = 0;
  unsigned int n_barcodes_unused = 0;
  size_t ctg_id, prevCtg_id = 0;

  if (! molecule_file.is_open()){
      std::cerr << "Error!  Cannot open file '" << Globals::molecule_file_name << "'.\n";
      exit(EXIT_FAILURE);
  }

  std::ofstream mapping_file;
  if (! Globals::mapping_file_name.empty()) {
    mapping_file.open(Globals::mapping_file_name);
  }
  long unsigned int n_lines;
  for (n_lines = 0; getline(molecule_file, molecule_line); ++n_lines){
    std::stringstream  splitstream(molecule_line);
    splitstream >> ctg >> beg_pos >> end_pos >> barcode >> nReads;
    auto pos = barcode_to_id.find(barcode);
    if (pos == barcode_to_id.end()) {
      barcode_id = barcode_to_id.size();
      barcode_to_id[barcode] = barcode_to_id.size();
    }
    else {
      barcode_id = pos->second;
    }
    if (ctg == prevCtg) {
      ctg_id = prevCtg_id;
    }
    else {
      auto pos = Globals::chrids.find(ctg);
      if (pos == Globals::chrids.end()) {
        std::cerr << "Error!  Contig '" << ctg << "' does not match any sequence name in the FASTA file.\n";
        exit(EXIT_FAILURE);
      }
      ctg_id     = pos->second;
      prevCtg    = ctg;
      prevCtg_id = ctg_id;
    }
    Interval molecule_interval (beg_pos, end_pos);
    for (ContigPart &contigPart: contigs[ctg_id].contigParts) {
      if (contigPart.get_overlap(molecule_interval) >= Globals::min_overlap) {
        // This is the condition to set a barcode to the beginning of a contig part
        if ((beg_pos >= contigPart.start) &&
            (beg_pos <= contigPart.start + std::min(Globals::window, contigPart.getSize() / 2)) &&
            (end_pos <= contigPart.start + contigPart.getSize() * (1 - Globals::beginning_ratio))) {
          contigPart.add_beg_molecule(barcode_id);
          ++n_barcodes_begin;
          if (! Globals::mapping_file_name.empty()) {
            mapping_file << "Barcode " << barcode << "\t" << ctg << "\t" << beg_pos << "\t" << end_pos << "\t" << nReads << " reads\t" << ctg << "\t" << contigPart.start << "\t" << contigPart.end << "\tbegin\n";
          }
        }
        // This is the condition to set a barcode to the end of a contig part
        else if ((end_pos <= contigPart.end) &&
            (beg_pos >= contigPart.start + contigPart.getSize() * Globals::beginning_ratio) && 
            (end_pos >= contigPart.end - std::min(Globals::window, contigPart.getSize() / 2))) {
          contigPart.add_end_molecule(barcode_id);
          ++n_barcodes_end;
          if (! Globals::mapping_file_name.empty()) {
            mapping_file << "Barcode " << barcode << "\t" << ctg << "\t" << beg_pos << "\t" << end_pos << "\t" << nReads << " reads\t" << ctg << "\t" << contigPart.start << "\t" << contigPart.end << "\tend\n";
          }
        }
        else {
          contigPart.add_other_molecule(barcode_id);
          ++n_barcodes_other;
          if (! Globals::mapping_file_name.empty()) {
            mapping_file << "Barcode " << barcode << "\t" << ctg << "\t" << beg_pos << "\t" << end_pos << "\t" << nReads << " reads\t" << ctg << "\t" << contigPart.start << "\t" << contigPart.end << "\tother\n";
          }
        }
      }
      else {
        ++n_barcodes_unused;
      }
    }
    if (n_lines % 1000000 == 0) std::cerr << n_lines << " lines read.\r" << std::flush;
  }
  std::cerr << n_lines << " lines read.\n";
  std::cerr << n_barcodes_begin << " anchored on the left part, " <<
    n_barcodes_end << " anchored on the right part, " <<
    (n_barcodes_begin + n_barcodes_end + n_barcodes_other) << " anchored in general, and " <<
    n_barcodes_unused << " unused.\n";
  std::cerr << barcode_to_id.size() << " different barcodes seen.\n";
  std::cerr << "Sorting barcodes:\n";
  for (size_t i = 0; i < contigs.size(); ++i) {
    Contig &contig = contigs[i];
    contig.sort_barcodes();
    if (i % 100 == 0) std::cerr << "\tContig #" << i << " / " << contigs.size() << ".\r" << std::flush;
  }
  std::cerr << "\tContig #" << contigs.size() << " / " << contigs.size() << ".\n";
}


Contig &get_contig (Contigs &contigs, NodeId &nodeId) {
  return contigs[nodeId.contigId];
}

ContigPart &get_contig_part (Contigs &contigs, NodeId &nodeId) {
  return contigs[nodeId.contigId].contigParts[nodeId.contigPartId];
}
