#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <map>

#include "fasta.h"


void complement (std::string &DNAseq) {
  reverse(DNAseq.begin(), DNAseq.end());
  for (size_t i = 0; i < DNAseq.length(); ++i){
    switch (DNAseq[i]) {
      case 'A':
        DNAseq[i] = 'T';
        break;    
      case 'C':
        DNAseq[i] = 'G';
        break;
      case 'G':
        DNAseq[i] = 'C';
        break;
      case 'T':
        DNAseq[i] = 'A';
        break;
      default:
        std::cerr << "Error! Unknown nucleotide '" << DNAseq[i] << "'.";
        exit(EXIT_FAILURE);
    }
  }
}

void read_fasta () {
  if (Globals::contigs_file_name.empty()) {
    std::cerr << "Error!  Input FASTA file is missing.\n";
    exit(EXIT_FAILURE);
  }
  std::ifstream ctg_fasta(Globals::contigs_file_name);
  if (! ctg_fasta.is_open()) {
    std::cerr << "Error!  Cannot open file '" << Globals::contigs_file_name << "'\n.";
    exit(EXIT_FAILURE);
  }
  std::string ctg_line, ctg_name, ctg_seq;
  while (getline(ctg_fasta, ctg_line)) {
    if (! ctg_line.empty()) {
      if (ctg_line[0] == '>') {
        if (! ctg_seq.empty()) {
          Globals::sequences[ctg_name] = ctg_seq;
          Globals::chrids[ctg_name] = Globals::chrs.size();
          Globals::chrs.push_back(ctg_name);
          Globals::chr_sizes.push_back(ctg_seq.size());
        }
        ctg_name = ctg_line.substr(1);
        ctg_seq.clear();
      }
      else {
        ctg_seq += ctg_line;
      }
    }
  }
  if (! ctg_seq.empty()) {
    Globals::sequences[ctg_name] = ctg_seq;
    Globals::chrids[ctg_name] = Globals::chrs.size();
    Globals::chrs.push_back(ctg_name);
    Globals::chr_sizes.push_back(ctg_seq.size());
  }
  ctg_fasta.close();
  std::cout << "Read " << Globals::chrs.size() << " contigs.\n";
  Globals::chrs.shrink_to_fit();
  Globals::chr_sizes.shrink_to_fit();
}

void write_fasta_sequence (int id, std::string &sequence, std::ofstream &file) {
  file << ">scaff_" << id << "\n";
  for (size_t i = 0; i <= sequence.length() / LINE_LEN; ++i) {
    file << sequence.substr(i * LINE_LEN, LINE_LEN) << "\n";
  }
}
