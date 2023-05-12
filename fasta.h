#ifndef FASTA_H
#define FASTA_H

#include <fstream>
#include <string>
#include <map>

const size_t line_len = 60;

void complement (std::string &DNAseq);
void read_fasta (std::map < std::string, std::string > &contigs, std::string &contig_fasta);
void write_fasta_sequence (int id, std::string &sequence, std::ofstream &file);

#endif
