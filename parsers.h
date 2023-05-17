#ifndef PARSERS_H
#define PARSERS_H

#include "molecule.h"
#include "interval.h"

void parse_contig_file();
void parse_molecule_file (Molecules &molecules);
void parse_contig_parts(RefIntervals &refIntervals);

#endif
