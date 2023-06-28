#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <cassert>
#include <iomanip>

#include "constants.h"
#include "contig.h"
#include "outliers_det.h"


// random sampling
// if seed = 0, seed is randomely set
void sampling (unsigned long n_sample, unsigned long max, std::vector < unsigned long > &id_sample) {
  std::vector < bool > used (max, false);
  unsigned long j = 0, nused = 0;
  unsigned long seed = Globals::seed;
  time_t t;
  if (seed == 0) seed = time(&t);
  // # samples is not less than sampling size
  if (n_sample >= max) {
    std::cerr << TAB << TAB << "No sampling.\n";
    for (unsigned long i = 0; i < n_sample; ++i) {
      id_sample[i] = i;
    }
    return;
  }
  std::cerr << TAB << "Sampling...\n";
  while (nused < n_sample) {
    srand(seed + j++);
    long id = rand() % max;
    if (! used[id]) {
      id_sample[nused] = id;
      used[id] = true;
      ++nused;
    }
  }
}


// input  -> X: data matrix (array), n: # of objects, d: # of dimensions, Xs:
// sample set, sid: sample indexes, ns: # of samples
// output -> result: an array of qsps
void qsp(std::vector < double > &X, unsigned long n, int d, unsigned long n_sample, std::vector < double > &score) {
  std::vector < unsigned long > id_sample (n_sample, 0);
  sampling(n_sample, n, id_sample);
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
    if (point % 10000 == 0) std::cerr << TAB << TAB << "Reading bin " << point << "/" << n << "\r" << std::flush;
  }
  std::cerr << TAB << TAB << "Reading bin " << n << "/" << n << "\n";
}

bool is_vector_constant (std::vector < double > &X) {
  double value = X.front();
  for (double x: X) {
    if (x != value) {
      return false;
	}
  }
  return true;
}

void normalize(std::vector < double > &X, unsigned long n) {
  double mean = 0.0;
  double sd   = 0.0;
  for (double x: X) {
    mean += x;
  }
  mean /= n;
  if (is_vector_constant(X)) {
    return;
  }
  for (double x: X) {
    sd += (x - mean) * (x - mean);
  }
  sd = sqrt(sd / (n - 1)); // unbiased SD
  for (unsigned long i = 0; i < n; i++) {
    X[i] /= sd;
  }
}


/*
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
      for (unsigned long i = 0; i < n; i++) {
        std::cerr << TAB << TAB << X[i * d + j] << "\n";
        X[i * d + j] = X[i * d + j] / sum;
	  }
    }
  }
}
*/

void compute_score (std::vector<double> &stat, unsigned long n, long d, unsigned long n_sample, std::vector<double> &score){
  score = std::vector < double > (stat.size(), 0);
  //normalize(stat, d, n);
  if (n_sample != 0) {
    n_sample = std::min(n_sample,n);
  }
  else {
    n_sample = n;
  }
  qsp(stat, n, d, n_sample, score);
}

double compute_threshold (std::vector < double > &scores, size_t n_elements) {
  double threshold = Globals::threshold;
  if (threshold != 0.0) {
    return threshold;
  }
  std::vector < double > sorted_scores = scores;
  std::sort(sorted_scores.begin(), sorted_scores.end());
  double q1 = sorted_scores[static_cast < unsigned long > (round(static_cast < double > (n_elements) * 0.25))];
  double q3 = sorted_scores[static_cast < unsigned long > (round(static_cast < double > (n_elements) * 0.75))];
  double iqr = q3 - q1;
  return q3 + 1.5 * iqr;
}

void detect_outliers (Molecule_stats &molecule_stats, std::vector < double > &scores, std::vector < bool > &outliers, size_t n_elements, size_t nchrs) {
  unsigned long cpt = 0;
  unsigned int d = 5;
  std::vector < std::vector < double > > values_tmp (d, std::vector < double > (n_elements));
  std::vector < double > all_values (d * n_elements);
  for (size_t chrid = 0; chrid < nchrs; ++chrid) {
    for (size_t i = 0; i < molecule_stats[chrid].size(); ++i) {
      values_tmp[0][cpt] = molecule_stats[chrid][i].coverage;
      values_tmp[1][cpt] = molecule_stats[chrid][i].length;
      values_tmp[2][cpt] = molecule_stats[chrid][i].read_density;
      values_tmp[3][cpt] = molecule_stats[chrid][i].start;
      values_tmp[4][cpt] = molecule_stats[chrid][i].end;
	  ++cpt;
    }
  }
  for (auto &X: values_tmp) {
    normalize(X, n_elements);
  }
  cpt = 0;
  for (unsigned int i = 0; i < n_elements; ++i) {
    for (unsigned int j = 0; j < d; ++j) {
      all_values[cpt++] = values_tmp[j][i];
    }
  }
  assert(cpt == d * n_elements);
  compute_score(all_values, n_elements, d, Globals::n_sample, scores);
  double threshold = compute_threshold(scores, n_elements);
  std::cerr << TAB << TAB << "Using a threshold of " << threshold << "\n";
  for (size_t i = 0; i < n_elements; ++i) {
    outliers[i] = (scores[i] <= threshold);
  }
}

void split (Molecule_stats &molecule_stats, Contigs &contigs, std::vector < double > &scores, std::vector < bool > &outliers, size_t n_elements, size_t nchrs) {
  std::ofstream output_file;
  std::ofstream scores_file;
  if (! Globals::scores_file_name.empty()) scores_file.open(Globals::scores_file_name, std::ofstream::out);
  if (! Globals::output_split_file_name.empty()) output_file.open(Globals::output_split_file_name, std::ofstream::out);
  unsigned int  n_bins_removed = 0;
  unsigned int  n_cuts         = 0;
  unsigned int  n_kept_chrs    = 0;
  unsigned long cpt            = 0;
  contigs.resize(Globals::chrs.size());
  for (size_t chrid = 0; chrid < nchrs; ++chrid) {
    unsigned int npos           = molecule_stats[chrid].size();
    bool         incut          = false;
    unsigned int bin_frag_start = 0;
    std::string &ctg            = Globals::chrs[chrid];
    for (unsigned int pos = 0; pos < npos; ++pos, ++cpt) {
      if (outliers[cpt]) {
        if (! incut) {
          // Check contig size
          int size = (pos - bin_frag_start) * Globals::window;
          if (size >= Globals::min_ctg_size) {
            // Do not print if we start with a cut
            if (pos > 0) {
              unsigned int start = bin_frag_start * Globals::window;
              unsigned int end   = pos * Globals::window;
              contigs[chrid].add_part(start, end);
			  ++n_cuts;
              // Output in BED format
              if (! Globals::output_split_file_name.empty()) {
                output_file << ctg << TAB << start << TAB << end << '\n';
              }
            }
          }
        }
        incut = true;
        ++n_bins_removed;
      }
      else if (incut) {
        incut          = false;
        bin_frag_start = pos;
      }
      if (! Globals::scores_file_name.empty()) scores_file << std::fixed << std::setprecision(6) <<
          ctg                         << TAB <<
          pos * Globals::window + 1   << TAB <<
          (pos + 1) * Globals::window << TAB <<
          scores[cpt] << '\n';
    }
    if (! incut) {
      unsigned int start = bin_frag_start * Globals::window;
      unsigned int end   = Globals::chr_sizes[chrid];
      contigs[chrid].add_part(start, end);
	  ++n_cuts;
      // Output in BED format
      if (! Globals::output_split_file_name.empty()) {
        output_file << ctg << TAB << start << TAB << end << '\n';
      }
    }
    if (! contigs[chrid].empty()) {
      ++n_kept_chrs;
	}
    if (chrid % 100 == 0) {
      std::cerr << TAB << "Splitting contig " << chrid << "/" << nchrs << "\r" << std::flush;
    }
  }
  assert(cpt == n_elements);
  std::cerr << TAB << "Splitting contig " << nchrs << "/" << nchrs << "\n";
  std::cerr << TAB << TAB << "Kept " << n_kept_chrs << " contigs, " << n_bins_removed << "/" << n_elements << " bins were removed, and " << n_cuts << " contig parts were created.\n";
  if (! Globals::output_split_file_name.empty()) output_file.close();
  if (! Globals::scores_file_name.empty()) scores_file.close();
}

void split_contigs (Molecule_stats &molecule_stats, Contigs &contigs) {
  size_t       nchrs          = Globals::chrs.size();
  size_t       n_elements     = 0;
  for (size_t chrid = 0; chrid < nchrs; ++chrid) {
    n_elements += molecule_stats[chrid].size();
  }
  std::vector < double > scores (n_elements);
  std::vector < bool > outliers (n_elements);
  detect_outliers(molecule_stats, scores, outliers, n_elements, nchrs);
  split(molecule_stats, contigs, scores, outliers, n_elements, nchrs);
}
