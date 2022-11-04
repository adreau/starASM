#ifndef DEF_INTERVAL
#define DEF_INTERVAL

struct Interval {
  
  unsigned long begin, end;

  Interval(unsigned long b, unsigned long e): begin(b), end(e) {}

  unsigned long int getSize () {
    return end - begin + 1;
  }
};

#endif
