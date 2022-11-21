#ifndef DEF_INTERVAL
#define DEF_INTERVAL

struct Interval {
  
  unsigned long begin, end;

  Interval (unsigned long b, unsigned long e): begin(b), end(e) {}

  unsigned long int getSize () {
    return end - begin + 1;
  }

  unsigned long int get_overlap (Interval &i) {
    int overlap = std::max(end, i.end) - std::min(begin, i.begin);
    if (overlap >= 0) return overlap;
    return 0;
  }

  unsigned long int get_distance (Interval &i) {
    return std::min(abs(begin - i.end), abs(end - i.begin));
  }

  void merge (Interval &i) {
    begin = std::min(begin, i.begin);
    end   = std::max(end,   i.end);
  }

  void unset () {
    begin = end = 0;
  }

  bool is_set () {
    return (begin != 0);
  }

  friend std::ostream &operator<< (std::ostream &out, Interval &e);
};

inline std::ostream &operator<< (std::ostream &out, Interval &i) {
  if (i.is_set()) {
    out << i.begin << "-" << i.end;
  }
  else {
    out << "--";
  }
  return out;
}

#endif
