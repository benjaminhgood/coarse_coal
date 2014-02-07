	#ifndef GENOME_HPP
#define GENOME_HPP

#include<vector>
#include "object_pool.hpp"

class Mutation{
    public:
        int location;
        //later on we will also have
        //double fitness effect
};

//inline bool operator<(Mutation const & m, const double p){ return m.location < p; }
//inline bool operator<(const double p, Mutation const & m){ return p < m.location; }
inline bool operator<(Mutation const & m1, Mutation const & m2){ return m1.location < m2.location; } 
//inline bool operator<=(Mutation const & m, const double p){ return m.location <= p; }
//inline bool operator<=(const double p, Mutation const & m){ return p <= m.location; }
inline bool operator<=(Mutation const & m1, Mutation const & m2){ return m1.location <= m2.location; }  
inline bool operator==(Mutation const & m1, Mutation const & m2){ return m1.location == m2.location; }
inline bool operator!=(Mutation const & m1, Mutation const & m2){ return m1.location != m2.location; }

typedef std::vector<Mutation> MutationList;
typedef SharedObjectPool<MutationList> GenomePool;

class Genome{
    public:

        int N;
        double s;
        double U;
        double R;
        double Un;
      
        double W;
        double NU;
        double NR;
        double NUn;

        int neutral_locus_location;   
        std::pair<int,int> marker_location;
        GenomePool genome_pool;

        Genome(int N, double s, double U, double R, double Un, double deltaLfrac): N(N), s(s), U(U), R(R), Un(Un) {

            W = exp(s);
            NU = N*U;
            NR = N*R;
            NUn = N*Un;

            int L = 1e+09;
            draw_site = create_random_int(0,L-1);


            neutral_locus_location = 0.5*L;
            marker_location.first = 0.5*(1-deltaLfrac)*L;
            marker_location.second = 0.5*(1+deltaLfrac)*L;

            genome_pool = GenomePool(2*N+2, MutationList());

            recombined = create_bernoulli(R);
            draw_num_mutations = create_poisson(U);
            draw_num_neutral_mutations = create_poisson(Un);

            draw_population_num_mutations = create_poisson(NU);
            draw_population_num_recombinants = create_poisson(NR);
            draw_population_num_neutral_mutations = create_poisson(NUn);

        };

        double draw_fitness_effect(Random & random){ return W; }; 
        decltype(create_random_int()) draw_site;
        
        decltype(create_bernoulli(R)) recombined;
        decltype(create_poisson(U)) draw_num_mutations;
        decltype(create_poisson(Un)) draw_num_neutral_mutations;

        decltype(create_poisson(NR)) draw_population_num_recombinants;
        decltype(create_poisson(NU)) draw_population_num_mutations;
        decltype(create_poisson(NUn)) draw_population_num_neutral_mutations;

        Mutation draw_mutation(Random & random){ 
            return Mutation{draw_site(random)}; 
        };
};

#endif
