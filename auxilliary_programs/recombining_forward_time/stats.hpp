#ifndef STATS_HPP
#define STATS_HPP

#include <random>
#include <functional>
#include <iostream>

typedef std::mt19937 Random;
static std::uniform_real_distribution<> sample_uniform(0,1);
static std::normal_distribution<> sample_normal();
static std::exponential_distribution<> sample_exponential();
static std::bernoulli_distribution sample_bernoulli(0.5); 

inline Random create_random(unsigned int seed=0, bool print_seed=false){
    if(seed==0){
        std::random_device rd;
        seed = rd();
    }
    if(print_seed){
        std::cout << seed << std::endl;
    }
    return Random(seed);
}

inline std::uniform_real_distribution<> create_uniform_real(double min, double max){ return std::uniform_real_distribution<>(min,max); }
inline std::uniform_int_distribution<> create_random_int(int min, int max){ return std::uniform_int_distribution<>(min,max);} 
inline std::uniform_int_distribution<> create_random_int(){ return std::uniform_int_distribution<>(); };
inline std::poisson_distribution<> create_poisson(double lam) { return std::poisson_distribution<>(lam);}
inline std::bernoulli_distribution create_bernoulli(double p){ return std::bernoulli_distribution(p);} 

inline int sample_poisson(Random & random, double mean){
    double L = exp(-mean);
    double p = sample_uniform(random);
    int k=0;
    while(p > L){
        ++k;
        p*=sample_uniform(random);
    }
    return k;
}

#endif
