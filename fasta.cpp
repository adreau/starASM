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

void read_fasta (std::map < std::string, std::string > &contigs, std::string &contig_fasta) {

    std::ifstream ctg_fasta(contig_fasta.c_str());
    if (! ctg_fasta.is_open()) {
        std::cerr << "Error!  Cannot open file '" << contig_fasta << "'" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string ctg_line, ctg_name, ctg_seq;

    while (getline(ctg_fasta, ctg_line)) {
        if (! ctg_line.empty()) {
            if (ctg_line[0] == '>') {
                if (! ctg_seq.empty()) {
                    contigs[ctg_name] = ctg_seq;
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
        contigs[ctg_name] = ctg_seq;
    }

    ctg_fasta.close();
    std::cout << "Read " << contigs.size() << " contigs.\n";
}

void write_fasta_sequence (int id, std::string &sequence, std::ofstream &file) {
    file << ">scaff_" << id << "\n";
    for (size_t i = 0; i <= sequence.length() / line_len; ++i) {
        file << sequence.substr(i * line_len, line_len) << "\n";
    }
}
