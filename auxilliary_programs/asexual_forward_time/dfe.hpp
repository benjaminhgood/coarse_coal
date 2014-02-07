#ifndef DFE_HPP
#define DFE_HPP

#include "stats.hpp"
#include "genome.hpp"

#include <vector>
#include <cmath>

class NeutralDFE {
    public: 
       NeutralDFE(double U): U(U) {};
       double get_fitness_effect(Random & random, Individual const & individual){ return 1;};
       void mutate_individual(Random & random, Individual const & individual) {};
    private:
       double U;
};

class DeltaDFE {
    public:
       DeltaDFE(double U, double s): U(U),s(s),W(std::exp(s)),L(std::exp(-U)){};

       double get_fitness_effect(Random & random, Individual const & individual) { 
           auto num_mutations = sample_precomputed_poisson(random,L);
           switch(num_mutations){
               case 0: return 1;
               case 1: return W;
               default: return std::pow(W,num_mutations);
           } };

       void mutate_individual(Random & random, Individual & individual) const {
           auto num_mutations = sample_precomputed_poisson(random,L);
           if(num_mutations > 0)
               individual.fitness *= std::pow(W,num_mutations);
       };

    private:
       double U;
       double s;
       double W;
       double L;        
};

class UniformDFE {
    public:
       UniformDFE(double U, double smin, double smax): W(std::exp(smax)), c(smin/smax), sample_num_mutations(U) {};

       void mutate_individual(Random & random, Individual & individual) {
           auto num_mutations = sample_num_mutations(random);
           double power = 0;
           for(int i=0;i<num_mutations;++i){
               power += c+(1-c)*sample_uniform(random);
           }
           individual.fitness *= std::pow(W,power);
       };

    private:
       std::poisson_distribution<> sample_num_mutations;
       double W;
       double c;
};


class TruncatedExponentialDFE {
    public:
       TruncatedExponentialDFE(double U, double smax, double c): W(std::exp(smax)), sample_num_mutations(U), p_thresh(1-std::exp(-c)) {};
    
       void mutate_individual(Random & random, Individual & individual) {
           auto num_mutations = sample_num_mutations(random);
           double power = 0;
           for(int i=0;i<num_mutations;++i){
               power += -log(1.0-sample_uniform(random)*p_thresh);
           }
           individual.fitness *= std::pow(W,power);
       };
    
    private:
        std::poisson_distribution<> sample_num_mutations;
        double W;
        double p_thresh;
 
};

class ExponentialDFE {
    public:
       ExponentialDFE(double U, double savg): W(std::exp(savg)), sample_num_mutations(U){};

       void mutate_individual(Random & random, Individual & individual) {
           auto num_mutations = sample_num_mutations(random);
           double power = 0;
           for(int i=0;i<num_mutations;++i){
               power += sample_exponential(random);
           }
           individual.fitness *= std::pow(W,power);
       };

    private:
       std::poisson_distribution<> sample_num_mutations;
       double W;
};

class GammaDFE{
    public:
        GammaDFE(double U, double savg, double alpha): W(std::exp(savg/alpha)), sample_num_mutations(U), sample_gamma(alpha, 1.0) {};
        
        void mutate_individual(Random & random, Individual & individual){
            auto num_mutations = sample_num_mutations(random);
            double power = 0;
            for(int i=0;i<num_mutations;++i){
                power += sample_gamma(random);
            }
            individual.fitness *= std::pow(W,power);
        };

    private:
        std::poisson_distribution<> sample_num_mutations;
        std::gamma_distribution<> sample_gamma;
        double W;
 };
        

class FiniteSitesDeltaDFE{
    public:
        FiniteSitesDeltaDFE(double U, double s, int L, double kL0);
        double get_fitness_effect(Random & random, Individual const & individual);
        void mutate_individual(Random & random, Individual & individual) {
            individual.fitness *= get_fitness_effect(random, individual);
        };
    private:
        double Utot;
        double s0;
        double kmu0;
        double mus0;
};

inline FiniteSitesDeltaDFE::FiniteSitesDeltaDFE(double U, double s, int L, double kL0) {
    Utot = U;
    double mu = Utot/L;
    s0 = s;
    mus0 = mu/s0;
    kmu0 = round(kL0*L)*mu;
}

inline double FiniteSitesDeltaDFE::get_fitness_effect(Random & random, Individual const & individual){
    double kmu = log(individual.fitness)*mus0+kmu0;
    return exp(s0*(sample_poisson(random, Utot-kmu)-sample_poisson(random,kmu)));
}

class SynonymousTrackedDFE{
    public:
       SynonymousTrackedDFE(double N, double U): U(U), mutated(create_bernoulli(U)), draw_num_mutations(create_poisson(N*U)), draw_mutated_index(create_uniform_int(0,((int) N)-1)) {};
       //bool mutated(Random & random) const { return (sample_uniform(random) < U); }; 
       decltype(create_bernoulli(0)) mutated;
       decltype(create_poisson(0)) draw_num_mutations;
       decltype(create_uniform_int(0,1)) draw_mutated_index; 
       Mutation get_mutation(Random & random, Individual const & individual) { return Mutation{label_generator.get_next_label(),0,0};};
    private:
       double U;
       LabelGenerator label_generator;       
};

#endif
