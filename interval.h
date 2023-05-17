#ifndef INTERVAL_H
#define INTERVAL_H

#include <vector>
#include <string>
#include <iostream>

#include "globals.h"

struct Interval {
  
  unsigned long start, end;

  Interval (unsigned long b, unsigned long e): start(b), end(e) {}

  unsigned long int getSize () {
    return end - start + 1;
  }

  unsigned long int get_overlap (Interval &i) {
    int overlap = std::max(end, i.end) - std::min(start, i.start);
    if (overlap >= 0) return overlap;
    return 0;
  }

  unsigned long int get_distance (const Interval &i) {
    return std::min(abs(start - i.end), abs(end - i.start));
  }

  void merge (Interval &i) {
    start = std::min(start, i.start);
    end   = std::max(end,   i.end);
  }

  void unset () {
    start = end = 0;
  }

  bool is_set () {
    return (start != 0);
  }

  friend std::ostream &operator<< (std::ostream &out, Interval &e);
};

inline std::ostream &operator<< (std::ostream &out, Interval &i) {
  if (i.is_set()) {
    out << i.start << "-" << i.end;
  }
  else {
    out << "--";
  }
  return out;
}


struct RefInterval: public Interval {

  std::string ref;
  bool strand;

  RefInterval (const std::string &r, bool s, unsigned long b, unsigned long e): Interval(b, e), ref(r), strand(s) {}

  bool can_merge (const RefInterval &next) {
    return ((ref == next.ref) && (strand == next.strand) && (get_distance(next) <= Globals::max_contig_distance));
  }

  friend std::ostream &operator<< (std::ostream &out, RefInterval &e);
};

inline std::ostream &operator<< (std::ostream &out, RefInterval &i) {
  if (i.is_set()) {
    out << (i.strand? '+': '-') << i.ref << " " << i.start << " " << i.end;
  }
  return out;
}

using RefIntervals = std::vector < RefInterval >;
using RefIntervalsSet = std::vector < std::vector < RefInterval > >;


void merge_close_contigs (RefIntervalsSet &refIntervalsSet);
void print_scaffold (RefIntervalsSet &refIntervalsSet);


#endif
