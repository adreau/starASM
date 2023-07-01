#ifndef READ_H
#define READ_H

#include <string>
#include "globals.h"


struct Read {
  unsigned long name;
  unsigned int chrid;
  unsigned long start, end;
  unsigned int mapq;
  Read (): name(0), chrid(0), start(0), end(0), mapq(0) {}
  Read (unsigned long n, unsigned int c, unsigned long s, unsigned long e, unsigned int m): name(n), chrid(c), start(s), end(e), mapq(m) {}
  bool is_solid () {
    return ((mapq >= Globals::min_mapq_solid) || ((end - start + 1) > Globals::min_len_solid));
  }
  unsigned long int get_distance (Read &i) {
    if (start   > i.end) return start - i.end;
    if (i.start > end)   return i.start - end;
    return 0;
  }
};

inline bool operator< (const Read &lhs, const Read &rhs) {
  if (lhs.chrid < rhs.chrid) {
    return true;
  }
  if (lhs.chrid > rhs.chrid) {
    return false;
  }
  return (lhs.start < rhs.start);
}

#endif
