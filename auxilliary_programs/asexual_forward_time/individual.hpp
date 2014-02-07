#ifndef INDIVIDUAL_HPP
#define INDIVIDUAL_HPP

#include <list>
#include <vector>
#include <memory>
#include <algorithm>
#include "stats.hpp"
#include "genome.hpp"

class Individual{
   public:
       double fitness;
       double location;
       std::shared_ptr<MutationList> tracked_mutations;
       Individual(){
           tracked_mutations = std::make_shared<MutationList>();
       };
       bool operator<(const double p) const { return location < p; };
       void add_mutation(Mutation m){
           if(tracked_mutations.unique()){
               tracked_mutations->push_back(m);
           }
           else{
               std::make_shared<MutationList>(*tracked_mutations).swap(tracked_mutations);
               tracked_mutations->push_back(m);
           }
       };
};

typedef std::vector<Individual> Population;

inline Population draw_random_sample(Random & random, Population const & population, int n){
    Population sample(n,population[0]);
    auto draw_random_index = create_random_int(0,population.size()-1);
    auto draw_random_individual = [&]()->Individual const &{ return population[draw_random_index(random)]; }; 
    std::generate(sample.begin(),sample.end(),draw_random_individual);
    return sample;
}
    


#endif
