#include "phylogeny.hpp"
#include <set>
#include <algorithm>

double calculate_pi(FrequencySpectrum const & f){
    // from Wakeley (2008)
    double pi=0;
    int n = f.size()+1;
    for(int i=1;i<n;++i){
        pi += f[i-1]*i*(n-i); 
    }
    pi /= (n*(n-1)/2);
    return pi;
}

double calculate_sn(FrequencySpectrum const & f){
    return f.sum();
}

FrequencySpectrum average_frequency_spectra(std::vector<FrequencySpectrum> const & frequency_spectra){
    if(frequency_spectra.size() < 1){
        return FrequencySpectrum(0.0,1);
    }
    else{
        FrequencySpectrum frequency_spectrum(0.0,frequency_spectra[0].size());
        for(auto & f : frequency_spectra){
            frequency_spectrum += f;
        }
        frequency_spectrum /= frequency_spectra.size();
        return frequency_spectrum;
    }
}

MarkerStatistics average_marker_statistics(std::vector<MarkerStatistics> const & marker_statistics_list){
    if(marker_statistics_list.size() < 1){
        return template_marker_statistics;
    }
    else{
        MarkerStatistics marker_statistics(template_marker_statistics);
        for(auto & m : marker_statistics_list){
            marker_statistics += m;
        }
        marker_statistics /= marker_statistics_list.size();
        return marker_statistics;
    }
}

MarkerStatistics calculate_marker_statistics(Population const & sample){
    MarkerStatistics marker_statistics(template_marker_statistics);
    for(auto & individual_1 : sample){
        for(auto & individual_2 : sample){
            double deltax = individual_1.marker_loci.first - individual_2.marker_loci.first;
            double deltay = individual_1.marker_loci.second - individual_2.marker_loci.second;
            marker_statistics[0] += deltax*deltax;
            marker_statistics[1] += deltay*deltay;
            marker_statistics[2] += deltax*deltax*deltax*deltax;
            marker_statistics[3] += deltay*deltay*deltay*deltay;
            marker_statistics[4] += deltax*deltax*deltay*deltay;
        }
    }
    marker_statistics /= sample.size()*(sample.size()-1);
    return marker_statistics;
}

FrequencySpectrum calculate_frequency_spectrum(Population const & sample){
    
    int n = sample.size();
    FrequencySpectrum frequency_spectrum(0.0,n-1);

    // merge polymorphic mutations into single list
    MutationList mutations;
    for(auto & individual : sample){
        mutations.insert(mutations.end(),individual.neutral_locus->begin(),individual.neutral_locus->end());
    } 

    // get unique mutations
    std::set<Mutation> unique_mutations(mutations.begin(),mutations.end());
    
    // construct frequency spectrum
    for(auto & unique_mutation : unique_mutations){
        int i = std::count(mutations.begin(),mutations.end(),unique_mutation);
        if(i<n) ++frequency_spectrum[i-1];
    }
    return frequency_spectrum;
}
