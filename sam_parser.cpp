#include <cassert>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include "htslib/htslib/sam.h"
#include "constants.h"
#include "globals.h"
#include "read.h"
#include "molecule.h"
#include "sam_parser.h"

bool is_mapped (unsigned int flag) {
  return ((flag & BAM_FUNMAP) == 0);
}

bool is_duplicate (unsigned int flag) {
  return ((flag & BAM_FDUP) != 0);
}

// Return the match length
unsigned int parse_cigar (std::string &cigar) {
  unsigned int match_len = 0;
  unsigned int number = 0;
  for (char c: cigar) {
    switch (c) {
      case '0':
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
        number = number * 10 + (c - '0');
        break;
      case 'M':
      case 'X':
      case '=':
      case 'D':
        match_len += number;
        number = 0;
        break;
      case 'I':
      case 'S':
      case 'P':
      case 'H':
        number = 0;
        break;
      default:
        assert(false);
    }
  }
  return match_len;
}

// Returns true iff:
//  - read is mapped
//  - flagged in proper pair
//  - mapq is higher than threshold
//  - barcode in BX tag is set
bool parse_main_line (std::string &line, unsigned int &chrid, unsigned long &start, unsigned long &end, unsigned int &mapq, std::string &barcode) {
  std::stringstream split_line(line);
  std::string value;
  const unsigned int no_value = -1;
  unsigned int flag = no_value;
  for (unsigned int i = 0; split_line >> value; ++i) {
    if (i == 1) {
      flag = std::stoi(value);
      if ((! is_mapped(flag)) || (is_duplicate(flag))) {
        return false;
      }
    }
    else if (i == 2) {
      if (value == "*") {
        return false;
      }
      assert(Globals::chrids.find(value) != Globals::chrids.end());
      chrid = Globals::chrids[value];
    }
    else if (i == 3) {
      start = std::stoul(value) - 1;
    }
    else if (i == 4) {
      mapq = std::stoi(value);
      if (mapq < Globals::min_mapq) {
        return false;
      }
    }
    else if (i == 5) {
    //hts_pos_t bam_endpos(const bam1_t *b);
      end = start + parse_cigar(value);
    }
    else if (i >= 11) {
      if (starts_with(value, BARCODE_FLAG)) {
        barcode = value.substr(6);
        //barcode = value.substr(BARCODE_FLAG.size());
      }
    }
  }
  if (barcode.empty()) {
    return false;
  }
  assert(chrid != no_value);
  assert(start != no_value);
  assert(end   != no_value);
  assert(mapq  != no_value);
  return true;
}

// Returns true iff:
//  - read is mapped
//  - flagged in proper pair
//  - mapq is higher than threshold
//  - barcode in BX tag is set
bool filter_line (bam1_t *aln, kstring_t *kbarcode) {
  unsigned int flag = aln->core.flag;
  if ((! is_mapped(flag)) || (is_duplicate(flag))) {
    return false;
  }
  if (aln->core.tid == -1) {
    return false;
  }
  if (aln->core.qual < Globals::min_mapq) {
    return false;
  }
  if (bam_aux_get_str(aln, BARCODE_FLAG, kbarcode) <= 0) {
    return false;
  }
  return true;
}

// Add the barcode count
// Returns true iff read passed the filters
bool add_barcode_count (Barcodes &barcodes, bam1_t *aln, kstring_t *kbarcode) {
  if (filter_line(aln, kbarcode)) {
    std::string barcode (ks_str(kbarcode), ks_len(kbarcode));
    ++barcodes.count_map[barcode];
    return true;
  }
  return false;
}

// Count the number of occurrences of each barcode
void count_barcodes (Barcodes &barcodes) {
  samFile          *fp           = hts_open(Globals::input_file_name.c_str(), "r");
  bam_hdr_t        *bam_hdr      = sam_hdr_read(fp);
  unsigned long int n_reads_kept = 0;
  unsigned long int cpt          = 0;
  std::cerr << "Reading BAM file for barcode counting...\n";
  bam1_t *aln = bam_init1();
  kstring_t kbarcode = KS_INITIALIZE;
  ks_initialize(&kbarcode);
  while (sam_read1(fp, bam_hdr, aln) > 0){
    if (add_barcode_count(barcodes, aln, &kbarcode)) {
      ++n_reads_kept;
    }
    ks_clear(&kbarcode);
    ++cpt;
    if (cpt % 10000000 == 0) {
      std::cerr << TAB << cpt << " lines read, " << n_reads_kept << " reads kept, using " << barcodes.count_map.size() << " barcodes.\r" << std::flush;
    }
  }
  ks_free(&kbarcode);
  bam_destroy1(aln);
  bam_hdr_destroy(bam_hdr);
  hts_close(fp);
  std::cerr << TAB << cpt << " lines read, " << n_reads_kept << " reads kept, using " << barcodes.count_map.size() << " barcodes.\n";
}

void add_barcode (Barcodes &barcodes, const std::string &name, unsigned int chrid, unsigned long start, unsigned long end, unsigned int mapq, const std::string &barcode) {
  auto it = barcodes.ids.find(barcode);
  if (it == barcodes.ids.end()) {
    return;
  }
  unsigned int id = it->second;
  assert(start < end);
  barcodes.reads[barcodes.current_offsets[id]] = Read(std::hash<std::string>{}(name), chrid, start, end, mapq);
  ++barcodes.current_offsets[id];
  assert(chrid < Globals::chrs.size());
}

void add_barcode_line (Barcodes &barcodes, bam1_t *aln, kstring_t *kbarcode) {
  if (filter_line(aln, kbarcode)) {
    std::string barcode (ks_str(kbarcode), ks_len(kbarcode));
    std::string name (bam_get_qname(aln), aln->core.l_qname);
    add_barcode(barcodes, name, aln->core.tid, aln->core.pos, bam_endpos(aln), aln->core.qual, barcode);
    ks_clear(kbarcode);
  }
}

void add_barcodes (Barcodes &barcodes) {
  samFile          *fp      = hts_open(Globals::input_file_name.c_str(), "r");
  bam_hdr_t        *bam_hdr = sam_hdr_read(fp);
  unsigned long int cpt     = 0;
  std::cerr << "Reading BAM file for barcode storing...\n";
  bam1_t *aln = bam_init1();
  kstring_t kbarcode = KS_INITIALIZE;
  ks_initialize(&kbarcode);
  while (sam_read1(fp, bam_hdr, aln) > 0){
    add_barcode_line(barcodes, aln, &kbarcode);
    ++cpt;
    if (cpt % 1000000 == 0) {
      std::cerr << TAB << cpt << " lines read.\r" << std::flush;
    }
  }
  barcodes.current_offsets.clear();
  ks_free(&kbarcode);
  bam_destroy1(aln);
  bam_hdr_destroy(bam_hdr);
  hts_close(fp);
  std::cerr << TAB << cpt << " lines read.\n";
}

void read_sam (Barcodes &barcodes) {
  if (Globals::input_file_name.empty()) {
    std::cerr << "Error!  Input BAM file is missing.\n";
    exit(EXIT_FAILURE);
  }
  count_barcodes(barcodes);
  barcodes.set_structure();
  add_barcodes(barcodes);
}
