#ifndef FASTA_H
#define FASTA_H

#include <fstream>
#include <string>
#include <map>

#include "constants.h"
#include "globals.h"

void complement (std::string &DNAseq);
void read_fasta ();
void write_fasta_sequence (int id, std::string &sequence, std::ofstream &file);

#endif
