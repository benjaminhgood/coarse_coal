#ifndef PHYLOGENY_HPP
#define PHYLOGENY_HPP

#include "stats.hpp"
#include "individual.hpp"
#include <valarray>
#include <list>

typedef std::valarray<double> FrequencySpectrum;

FrequencySpectrum calculate_frequency_spectrum(Population const & sample);
double calculate_pi(Population const & sample);
double calculate_pi(Individual const & individual_1, Individual const & individual_2);
double calculate_pi(FrequencySpectrum const & frequency_spectrum);
double calculate_sn(FrequencySpectrum const & frequency_spectrum);
std::list<double> calculate_sparse_pis(Population const & sample);

template <class Container>
inline typename Container::value_type mean(Container const & c) 
{
    typename Container::value_type init(c.front());
    init-=init;
    for(auto & item : c){
        init+=item;
    }
    init /= c.size();
    return init;
}


template <class Container>
inline typename Container::value_type mean_squared(Container const & c) 
{
    typename Container::value_type init(c.front());
    init-=init;
    for(auto & item : c){
        auto temp(item);
        temp*=item;
        init+=temp;
    }
    init /= c.size();
    return init;
}

#endif
