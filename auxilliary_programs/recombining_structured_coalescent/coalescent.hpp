#ifndef COALESCENT_HPP
#define COALESCENT_HPP

#include "stats.hpp"
#include <valarray>
typedef std::valarray<double> FrequencySpectrum;
FrequencySpectrum simulate_frequency_spectrum(Random & random, double Ns, double NU, double NR, int n, int total_runs);

#endif 
