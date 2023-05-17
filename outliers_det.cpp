#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <cassert>
#include <iomanip>

#include "constants.h"
#include "outliers_det.h"


// random sampling
// if seed = 0, seed is randomely set
void sampling (unsigned long n_sample, unsigned long max, long seed, std::vector < unsigned long > &id_sample) {
  std::vector < bool > used (max, false);
  unsigned long j = 0, nused = 0;
  time_t t;

  // # samples is not less than sampling size
  if (n_sample >= max) {
    std::cerr << "\tNo sampling.\n";
    for (unsigned long i = 0; i < n_sample; ++i) {
      id_sample[i] = i;
    }
    return;
  }

  std::cerr << "\tSampling...\n";
  if (seed == 0)
    seed = time(&t);

  while (nused < n_sample) {
    srand(seed + j++);
    long id = rand() % max;
    //long id = rand() % (max + 1);
    if (! used[id]) {
      id_sample[nused] = id;
      used[id] = true;
      ++nused;
    }
  }
  std::cerr << "\t... done.\n";
}


// input  -> X: data matrix (array), n: # of objects, d: # of dimensions, Xs:
// sample set, sid: sample indexes, ns: # of samples
// output -> result: an array of qsps
void qsp(std::vector < double > &X, unsigned long n, int d, unsigned long n_sample, long seed, std::vector < double > &score) {
  std::vector < unsigned long > id_sample (n_sample, 0);

  sampling(n_sample, n, seed, id_sample);

  // compute the outlierness score qsp for each data point
  for (unsigned long point = 0; point < n; point++) {
    double res = 0;
    for (long unsigned i = 0; i < n_sample; i++) {
      if (point != id_sample[i]) {
        double sum = 0;
        for (int j = 0; j < d; j++) {
          sum += fabs(X[point * d + j] - X[id_sample[i] * d + j]) *
                 fabs(X[point * d + j] - X[id_sample[i] * d + j]);
        }
        if (res == 0)
          res = sum;
        else if (sum < res && sum > 0)
          res = sum;
      }
    }
    score[point] = sqrt(res);
    if (point % 10000 == 0) std::cerr << "\t" << point << "/" << n << "\r" << std::flush;
  }
  std::cerr << "\t" << n << "/" << n << "\n";
}


// normalization of data (divide by SDs for each dimension)
void normalize(std::vector < double > &X, unsigned long n, int d) {
  std::vector < double > X_means (d, 0);

  for (int j = 0; j < d; j++) {
    double sum = 0;
    for (unsigned long int i = 0; i < n; i++) {
      sum += X[i * d + j];
    }
    X_means[j] = sum / n;
  }

  for (int j = 0; j < d; j++) {
    bool flag = false;
    double sum = X[j];
    for (unsigned long int i = 1; i < n; i++) {
      if (sum != X[i * d + j]) {
        flag = true;
        break;
      }
    }
    if (flag) {
      sum = 0;
      for (unsigned long i = 0; i < n; i++)
        sum += (X[i * d + j] - X_means[j]) * (X[i * d + j] - X_means[j]);
      sum = sqrt(sum / (n - 1)); // unbiased SD
      for (unsigned long i = 0; i < n; i++)
        X[i * d + j] = X[i * d + j] / sum;
    }
  }
}

void compute_score(std::vector<double> &stat, unsigned long n, long d, unsigned long n_sample, long seed, std::vector<double> &score){

  score = std::vector < double > (stat.size(), 0);
  normalize(stat, d, n);

  if(n_sample != 0) {
    n_sample = std::min(n_sample,n);
  }else{
    n_sample = n;
  }
  qsp(stat, n, d, n_sample, seed, score);
}


void detect_outliers(Molecule_stats &molecule_stats) {

  std::ofstream output_file;
  std::ofstream scores_file;
  if (! Globals::scores_file_name.empty()) scores_file.open(Globals::scores_file_name, std::ofstream::out);
  if (! Globals::output_split_file_name.empty()) output_file.open(Globals::output_split_file_name, std::ofstream::out);

  size_t nchrs = Globals::chrs.size();
  size_t n_elements = 0;
  int seed = 0;
  for (size_t chrid = 0; chrid < nchrs; ++chrid) {
    n_elements += molecule_stats[chrid].size();
  }
  int d = 5;
  std::vector < double > all_values (d * n_elements);
  std::vector < double > score_all_values (n_elements);
  size_t cpt = 0;
  for (size_t chrid = 0; chrid < nchrs; ++chrid) {
    for (size_t i = 0; i < molecule_stats[chrid].size(); ++i) {
      all_values[cpt++] = molecule_stats[chrid][i].coverage;
      all_values[cpt++] = molecule_stats[chrid][i].length;
      all_values[cpt++] = molecule_stats[chrid][i].read_density;
      all_values[cpt++] = molecule_stats[chrid][i].start;
      all_values[cpt++] = molecule_stats[chrid][i].end;
    }
  }
  assert(cpt == d * n_elements);

  compute_score(all_values, n_elements, d, Globals::n_sample, seed, score_all_values);

  cpt = 0;
  for (size_t chrid = 0; chrid < nchrs; ++chrid) {
    unsigned int npos           = molecule_stats[chrid].size();
    bool         incut          = false;
    unsigned int bin_frag_start = 0;
    std::string &ctg    = Globals::chrs[chrid];

    for (unsigned int pos = 0; pos < npos; ++pos, ++cpt) {
      if (score_all_values[cpt] > Globals::threshold) {
        if (! incut) {
          // Check contig size
          int size = (pos - bin_frag_start) * Globals::window;
          if (size >= Globals::min_ctg_size) {
            // Do not print if we start with a cut
            // Output in BED format
            if (pos > 0) {
              if (! Globals::output_split_file_name.empty()) output_file << ctg << '\t' << bin_frag_start * Globals::window << '\t' << pos * Globals::window << '\n';
            }
          }
        }
        incut = true;
      }
      else if (incut) {
        incut          = false;
        bin_frag_start = pos;
      }
      if (! Globals::scores_file_name.empty()) scores_file << std::fixed << std::setprecision(6) <<
          ctg                         << tab <<
          pos * Globals::window + 1   << tab <<
          (pos + 1) * Globals::window << tab <<
          score_all_values[cpt] << '\n';
    }
    std::cerr << "Analyzing contig #" << chrid << "/" << nchrs << ".\r" << std::flush;
    if (! incut) {
      // Output in BED format
      if (! Globals::output_split_file_name.empty()) output_file << ctg << tab << bin_frag_start * Globals::window << tab << Globals::chr_sizes[chrid] << '\n';
    }
  }
  std::cerr << "Analyzing contig #" << nchrs << "/" << nchrs << ".\n";

  if (! Globals::output_split_file_name.empty()) output_file.close();

  if (! Globals::scores_file_name.empty()) scores_file.close();
}
