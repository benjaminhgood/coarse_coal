#ifndef EVOLVE_HPP
#define EVOLVE_HPP

#include "stats.hpp"
#include "individual.hpp"
#include "dfe.hpp"
#include "phylogeny.hpp"

class Results {
    public:
        MutationList fixed_mutations;
        std::vector<double> avg_fitnesses;
        std::vector<int> record_times;
        std::vector<FrequencySpectrum> frequency_spectra;
        std::list<double> pis;
};


template <class UntrackedDFE, class TrackedDFE>
inline Results evolve_population(Random & random, double floatN, UntrackedDFE & dfe, TrackedDFE & tracked_dfe, int total_samples, int n){

    Results results;
    int N = floatN;    
    Population population(N,Individual());
    Population new_population(population);    

    // initialize population
    double total_fitness = 0;
    for(auto & individual : population){
        individual.fitness = 1.0;
        individual.add_mutation(tracked_dfe.get_mutation(random));
        total_fitness += individual.fitness;
        individual.location = total_fitness;
    }

    results.avg_fitnesses.push_back(total_fitness/N);
    results.record_times.push_back(0);

    bool in_equilibrium = false;
    int deltat = 0;
    // loop over desired number of generations to evolve
    for(int t=0,num_samples=0;num_samples < total_samples;++t){
        
        // create a functor for drawing parents
        auto draw_parent = [&population,total_fitness](Random & random)->Individual const & { return *std::lower_bound(population.begin(),population.end(),total_fitness*sample_uniform(random)); };
        
        // construct the next generation 
        total_fitness = 0;
        for(auto & individual : new_population){
            individual=draw_parent(random); 
            //individual.fitness *= dfe.get_fitness_effect(random, individual);
            dfe.mutate_individual(random,individual);
            /*if(tracked_dfe.mutated(random)){
                individual.add_mutation(tracked_dfe.get_mutation(random));
            }*/
            total_fitness+=individual.fitness;
            individual.location = total_fitness;
        }
        std::swap(population,new_population);

        // add tracked mutations
        for(int i=0, num_tracked_mutations=tracked_dfe.draw_num_mutations(random);i<num_tracked_mutations;++i){
            Individual & mutated_individual = population[tracked_dfe.draw_mutated_index(random)];
            Mutation new_mutation = tracked_dfe.get_mutation(random, mutated_individual);
            new_mutation.W_background = mutated_individual.fitness / avg_fitness;
            mutated_individual.add_mutation(new_mutation);
        }
         
        // record stuff if is the right time
        if(in_equilibrium && (t % deltat == 0)){
            ++num_samples;
            results.record_times.push_back(t);
            results.avg_fitnesses.push_back(total_fitness/N);
            for(int i=0;i<100;++i){ 
                auto sample = draw_random_sample(random,population,n);
                results.frequency_spectra.push_back(calculate_frequency_spectrum(sample));
                //results.pis.splice(results.pis.end(),calculate_sparse_pis(sample));
                results.pis.push_back(calculate_pi(sample));
            }
        }
        
        // check for a fixed mutation!
        if(!population.front().tracked_mutations->empty()){
            // can't have a fixed mutation without a mutation in the first lineage
            auto first_mutation = population.front().tracked_mutations->front();
            // create a functor for testing whether a lineage shares that mutation
            auto same_first_mutation = [first_mutation](Individual const & individual)
            { 
               return !individual.tracked_mutations->empty() && individual.tracked_mutations->front() == first_mutation; 
            };
            if(std::all_of(population.begin(),population.end(),same_first_mutation)){
                // if all lineages share that first mutation,
                // record it
                results.fixed_mutations.push_back(first_mutation); 
                // and remove it from the population
                for(auto & individual : population){
                    if(same_first_mutation(individual)){
                        individual.tracked_mutations->erase(individual.tracked_mutations->begin());
                    } 
                }
                // we've decided that we are in equilibrium if a mutation fixes 
                if(!in_equilibrium){
                    in_equilibrium = true;
                    deltat = t;
                    results.record_times.push_back(t);
                    results.avg_fitnesses.push_back(total_fitness/N);
                    //std::cout << "Reached equilibrium " << t << " " << total_fitness/N << std::endl;
                }
            } 
        }
    }
    //std::cout << "Finished " << results.record_times.back() << " " << results.avg_fitnesses.back() << std::endl;
    return results;
}

#endif 
