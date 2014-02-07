#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <math.h>
#include "stats.hpp"
#include "genome.hpp"
#include "individual.hpp"
#include "phylogeny.hpp"

class Measurement {
    public:
        int time;
        double avg_fitness;
        MarkerStatistics marker_statistics;
        FrequencySpectrum frequency_spectrum;
};

class Results {
    public:
        int equilibrium_time;
        double equilibrium_avg_fitness;
        MutationList fixed_mutations;
        std::vector<Measurement> measurements;
};


Results evolve_population(Random & random, double floatN, Genome & genome, int t_max, int n);
MutationList remove_fixed_genomic_mutations(Population & population);

int main(int argc, char * argv[]){
    if(argc < 9){
        std::cout << "usage: " << argv[0] << " t_max n N s U R Un deltaR" << std::endl;
        return 1;
    }
    else{
        unsigned long t_max = atol(argv[1]);
        unsigned long n = atol(argv[2]);
        double N = atof(argv[3]);
        double s = atof(argv[4]);
        double U = atof(argv[5]);
        double R = atof(argv[6]);
        double Un = atof(argv[7]);
        double deltaLfrac = atof(argv[8]);

        //Random random = create_random(42);
        Random random = create_random();
        Genome genome(N,s,U,R,Un,deltaLfrac);
        Results results = evolve_population(random,N,genome,t_max,n);
        
        if(results.measurements.empty()){
            std::cout << "Did not evolve long enough to make a measurement!" << std::endl;
        }
        else{
            auto avg_frequency_spectrum = FrequencySpectrum(0.0, results.measurements.back().frequency_spectrum.size());
            auto avg_marker_statistics = template_marker_statistics;

            for(auto & measurement : results.measurements){
                avg_frequency_spectrum += measurement.frequency_spectrum;
                avg_marker_statistics += measurement.marker_statistics;                
            }
            avg_frequency_spectrum /= results.measurements.size();
            avg_marker_statistics /= results.measurements.size();

            double pi = calculate_pi(avg_frequency_spectrum);
            double Sn = calculate_sn(avg_frequency_spectrum);
        
            std::cout << log(results.measurements.back().avg_fitness/results.equilibrium_avg_fitness)/(results.measurements.back().time-results.equilibrium_time) << " " << pi << " " << Sn << std::endl;
            std::cout << (avg_marker_statistics[0]+avg_marker_statistics[1])/2 << " ";
            std::cout << (avg_marker_statistics[2]+avg_marker_statistics[3])/2 << " ";
            std::cout << avg_marker_statistics[4] << std::endl;
            std::cout << results.fixed_mutations.size() << std::endl; 

            for(double fi : avg_frequency_spectrum){
                std::cout << fi << " ";
            }
            std::cout << std::endl; 
            return 0;
        }
    }
}

Results evolve_population(Random & random, double floatN, Genome & genome, int t_max, int n){

    Results results;
    int N = floatN;
    int total_samples = t_max/N;
    auto draw_index = create_random_int(0,N-1);

    Population population(N, Individual(genome)); // the current population of individuals
    Population new_population(N, population.front()); // used for constructing the new generation     
    Population sample(n, population.front()); // for sampling diversity  
    
    // initialize population
    double total_fitness = 0;
    for(auto & individual : population){
        individual.mutate_neutral_locus(random, genome);
        individual.fitness = 1.0;
        total_fitness += individual.fitness;
        individual.location = total_fitness;
    }

    bool in_equilibrium = false;
    int deltat = 0;

    for(int t=0;t< t_max;++t){
        
        auto draw_parent = [&,total_fitness](Random & r)->Individual & { return *std::lower_bound(population.begin(),population.end(),total_fitness*sample_uniform(r)); };
        
        for(auto & individual : new_population){
            individual = draw_parent(random);
            individual.mutate_marker_loci(random, genome);
        }

        // do recombination 
        for(int i=0,num_recombinants = genome.draw_population_num_recombinants(random);i<num_recombinants;++i){
            new_population[draw_index(random)].recombine(random, genome, draw_parent(random)); 
        }

        // add genomic mutations
        for(int i=0,num_mutations = genome.draw_population_num_mutations(random);i<num_mutations;++i){
            new_population[draw_index(random)].mutate_coding_region(random, genome); 
        }

        // add neutral loci mutations
        for(int i=0,num_mutations = genome.draw_population_num_neutral_mutations(random);i<num_mutations;++i){
            new_population[draw_index(random)].mutate_neutral_locus(random, genome); 
        }

        total_fitness = 0;        
        for(auto & individual : new_population){
            total_fitness += individual.fitness;
            individual.location = total_fitness;
        }
        std::swap(population,new_population);

 
        if(t % 100 == 0){
            auto new_fixed_mutations = remove_fixed_genomic_mutations(population);
            results.fixed_mutations.insert(results.fixed_mutations.end(), new_fixed_mutations.begin(), new_fixed_mutations.end());
        }
        // record stuff if necessary
        if(in_equilibrium && (t % deltat == 0)){
            for(int i=0;i<100;i++){
                // draw sample
                for(auto & individual : sample){
                    individual = population[draw_index(random)];
                }
                results.measurements.push_back(Measurement{t, total_fitness/floatN*std::pow(genome.W,results.fixed_mutations.size()), calculate_marker_statistics(sample), calculate_frequency_spectrum(sample)});
            }
        }


        // check to see if a mutation has fixed at the neutral locus
        if(!population.front().neutral_locus->empty() ){ 
            // a mutation can only fix if the first individual has it
            auto first_mutation = population.front().neutral_locus->front();
            auto same_first_mutation = [first_mutation](Individual const & individual){return (!individual.neutral_locus->empty()) && (individual.neutral_locus->front() == first_mutation);};
 
            if(std::all_of(++population.begin(),population.end(),same_first_mutation)){
                for(auto & individual : population){
                    if(same_first_mutation(individual)) 
                        individual.neutral_locus->erase(individual.neutral_locus->begin());
                }
                //results.fixed_mutations.push_back(first_mutation);
                if(!in_equilibrium){
                    in_equilibrium = true;
                    results.equilibrium_time = t;
                    results.equilibrium_avg_fitness = total_fitness/floatN*std::pow(genome.W, results.fixed_mutations.size());
                    deltat = t;
                }
            }
        }


    }
    return results;
}


MutationList remove_fixed_genomic_mutations(Population & population){
    MutationList fixed_mutations(population.front().mutations->begin(),population.front().mutations->end());
    MutationList new_fixed_mutations = fixed_mutations;

    //std::cout << fixed_mutations.size() << std::endl;

    for(auto & individual : population){
        auto new_end = std::set_intersection(fixed_mutations.begin(),fixed_mutations.end(), individual.mutations->begin(),individual.mutations->end(), new_fixed_mutations.begin());
        new_fixed_mutations.resize(new_end-new_fixed_mutations.begin());
        std::swap(fixed_mutations, new_fixed_mutations);
        if(fixed_mutations.empty()) break;
    }
    if(!fixed_mutations.empty()){
        //std::cout << "Removing genomic mutations!" << std::endl;
            for(auto & mutation : fixed_mutations){
                for(auto & individual : population){
                    auto mutation_ptr = std::lower_bound(individual.mutations->begin(), individual.mutations->end(), mutation);
                    if(mutation_ptr != individual.mutations->end() && *mutation_ptr == mutation){
                        individual.mutations->erase(mutation_ptr);
                    }
                }
            }
    }
    return fixed_mutations;
}







