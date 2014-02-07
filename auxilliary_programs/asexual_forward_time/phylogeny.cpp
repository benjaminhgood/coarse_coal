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

FrequencySpectrum calculate_frequency_spectrum(Population const & sample){
    
    int n = sample.size();
    FrequencySpectrum frequency_spectrum(0.0,n-1);

    // merge polymorphic mutations into single list
    MutationList mutations;
    for(auto & individual : sample){
        mutations.insert(mutations.end(),individual.tracked_mutations->begin(),individual.tracked_mutations->end());
    } 

    std::set<Mutation> unique_mutations(mutations.begin(),mutations.end());
    std::multiset<Mutation> sorted_mutations(mutations.begin(),mutations.end());
    
    // construct frequency spectrum
    for(auto & unique_mutation : unique_mutations){
        auto i = sorted_mutations.count(unique_mutation); 
        if(i<n) ++frequency_spectrum[i-1];
    }
    return frequency_spectrum;
}

double calculate_pi(Population const & sample){
    auto n = sample.size();
    double pi = 0;
    for(int i=0; i<n; ++i){
        for(int j=i+1; j<n; ++j){
            pi += calculate_pi(sample[i],sample[j]);
        }
    }
    pi = pi/(n*(n-1)/2);
    return pi;
}

std::list<double> calculate_sparse_pis(Population const & sample){
    std::list<double> pis;
    for(int i=0,j=1,n=sample.size();j<n;i+=2,j+=2){
        pis.push_back(calculate_pi(sample[i],sample[j]));
    }
    return pis;
}

double calculate_pi(Individual const & individual_1, Individual const & individual_2){
    // merge polymorphic mutations into a single list
    MutationList mutations;
    mutations.insert(mutations.end(),individual_1.tracked_mutations->begin(),individual_1.tracked_mutations->end());
    mutations.insert(mutations.end(),individual_2.tracked_mutations->begin(),individual_2.tracked_mutations->end());
    // get unique mutations
    std::set<Mutation> unique_mutations(mutations.begin(),mutations.end());
    
    return 2*unique_mutations.size()-mutations.size();
}

