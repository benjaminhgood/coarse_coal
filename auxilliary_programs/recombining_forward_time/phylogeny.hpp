#ifndef PHYLOGENY_HPP
#define PHYLOGENY_HPP

#include "stats.hpp"
#include "individual.hpp"
#include <valarray>

typedef std::valarray<double> FrequencySpectrum;
typedef std::valarray<double> MarkerStatistics;
static MarkerStatistics template_marker_statistics(0.0,5); 
// a[0]=deltax, a[1]=deltay, a[2]=deltax^2, a[3]=deltay^2, a[4]=deltax*deltay

FrequencySpectrum calculate_frequency_spectrum(Population const & population);
MarkerStatistics calculate_marker_statistics(Population const & population);

FrequencySpectrum average_frequency_spectra(std::vector<FrequencySpectrum> const & frequency_spectra);
MarkerStatistics average_marker_statistics(std::vector<MarkerStatistics> const & marker_statistics);

double calculate_pi(FrequencySpectrum const & frequency_spectrum);
double calculate_sn(FrequencySpectrum const & frequency_spectrum);

#endif
