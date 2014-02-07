#ifndef INDIVIDUAL_HPP
#define INDIVIDUAL_HPP

#include <vector>
#include <memory>
#include <algorithm>
#include "stats.hpp"
#include "genome.hpp"

class Individual{
   public:
       double location; // used for intrusive find while drawing parents
       double fitness; // used for memoization

       decltype(GenomePool().allocate()) mutations;       
       std::shared_ptr<MutationList> neutral_locus;
       std::pair<double,double> marker_loci;

       Individual() {};
       Individual(Genome & reference_genome): mutations(reference_genome.genome_pool.allocate()), neutral_locus(std::make_shared<MutationList>()), marker_loci{0,0} {};       

       void mutate_coding_region(Random & r, Genome & reference_genome);
       void mutate_marker_loci(Random & r, Genome & reference_genome);
       void mutate_neutral_locus(Random & r, Genome & reference_genome);

       void recombine(Random & random, Genome & reference_genome, Individual & other_parent);

       void recalculate_fitness(Genome & reference_genome){ fitness = pow(reference_genome.W, mutations->size()); };

       bool operator<(const double p) const { return location < p; };        
};

inline void Individual::mutate_coding_region(Random & random, Genome & reference_genome){
    auto new_mutation = reference_genome.draw_mutation(random);
    auto new_mutations = reference_genome.genome_pool.allocate();
    auto insertion_point = std::upper_bound(mutations->begin(),mutations->end(),new_mutation);
    new_mutations->assign(mutations->begin(), insertion_point);
    new_mutations->push_back(new_mutation);
    if(mutations->end() != insertion_point)
        new_mutations->insert(new_mutations->end(), insertion_point, mutations->end());

    mutations.swap(new_mutations);
    fitness *= reference_genome.W;
}

inline void Individual::mutate_neutral_locus(Random & random, Genome & reference_genome){
    std::make_shared<MutationList>(*neutral_locus).swap(neutral_locus);
    neutral_locus->push_back(reference_genome.draw_mutation(random));   
}

inline void Individual::mutate_marker_loci(Random & random, Genome & reference_genome){
    sample_bernoulli(random) ? marker_loci.first+=1 : marker_loci.first-=1;
    sample_bernoulli(random) ? marker_loci.second+=1 : marker_loci.second-=1;
}

inline void Individual::recombine(Random & random, Genome & reference_genome, Individual & other_parent){
    // do coding region
    auto new_mutations = reference_genome.genome_pool.allocate();
    
    int breakpoint = reference_genome.draw_site(random);
    auto break_point = std::upper_bound(mutations->begin(),mutations->end(), Mutation{breakpoint});
    new_mutations->assign(mutations->begin(), break_point);
    break_point = std::upper_bound(other_parent.mutations->begin(),other_parent.mutations->end(), Mutation{breakpoint});
    if( break_point != other_parent.mutations->end() )
        new_mutations->insert(new_mutations->end(), break_point, other_parent.mutations->end());
    mutations.swap(new_mutations);

    recalculate_fitness(reference_genome);

    // now do marker loci
    if(breakpoint < reference_genome.marker_location.second){
        if(breakpoint < reference_genome.marker_location.first){
                marker_loci = other_parent.marker_loci;
        }
        else{
            marker_loci.second = other_parent.marker_loci.second;    
        }
    }

    // do neutral locus
    if(breakpoint < reference_genome.neutral_locus_location)
        neutral_locus = other_parent.neutral_locus;
}

typedef std::vector<Individual> Population;
    
inline Population draw_random_sample(Random & random, Population const & population, int n){
    Population sample(n,population[0]);
    auto draw_random_index = create_random_int(0,population.size()-1);
    auto draw_random_individual = [&]()->Individual const &{ return population[draw_random_index(random)]; }; 
    std::generate(sample.begin(),sample.end(),draw_random_individual);
    return sample;
}

#endif
