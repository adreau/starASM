#ifndef DEF_INTERVAL
#define DEF_INTERVAL

struct Interval {
  
  unsigned long begin, end;

  Interval (unsigned long b, unsigned long e): begin(b), end(e) {}

  unsigned long int getSize () {
    return end - begin + 1;
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
    return (begin == 0);
  }
};

#endif
