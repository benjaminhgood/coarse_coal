#include <iostream>
#include <vector>
#include <algorithm>
#include <random>

// For random number generation
typedef std::mt19937 Random;
typedef std::poisson_distribution<> poisson_distribution;
typedef std::uniform_int_distribution<> uniform_int_distribution;
std::uniform_real_distribution<> sample_uniform(0,1);
std::exponential_distribution<> sample_exponential(1);

class FitnessClass{ // Tracks the ancestral lineages that currently reside in this class
    public:
        double coalescence_rate; // rate at which two lineages in this class coalesce
        double mutation_rate_x; // rate at which lineages leave this class due to a mutation event
        double mutation_rate_y; // rate at which lineages leave this class due to a mutation event        
        std::vector<int> lineages; // number of present day descendents of each lineage

        FitnessClass(double coalescence_rate=0, double mutation_rate_x=0, double mutation_rate_y=0, int size_hint = 1): coalescence_rate(coalescence_rate), mutation_rate_x(mutation_rate_x), mutation_rate_y(mutation_rate_y) { lineages.reserve(size_hint); }
};

typedef std::vector<std::vector<FitnessClass>> FitnessDistribution;

class Event{
    public:
        enum EventType { MUTATION_X, MUTATION_Y, COALESCENCE};
        double weight; // used for intrusive search
        EventType type;
        int kx;
        int ky; // the fitness class where this event occurs
        bool operator<(const double p) const { return weight < p; } // used for intrusive search
};

// The primary function in this file. Simulates the synonymous site frequency spectrum.
std::vector<double> simulate_frequency_spectrum(Random & random, double Ns1, double NUd1, double Ns2, double NUd2, int sample_size, int num_samples);


// Initializes the population fitness distribution according to the supplied parameters
// Doesn't create all the fitness classes initially, but only up to 4*lambda+1
FitnessDistribution create_fitness_distribution(double Ns1, double lambda1, double Ns2, double lambda2, int size_hint);

// Extends the fitness distribution to the new k_max if necessary
void extend_fitness_distribution_x(FitnessDistribution & fitness_distribution, int k_max);
void extend_fitness_distribution_y(FitnessDistribution & fitness_distribution, int k_max);


// Main function. Reads in parameters, calls simulate_frequency_spectrum, and prints results.
int main(int argc, char* argv[]){
    if(argc < 7){
        std::cout << "usage: " << argv[0] << " Ns1 NUd1 Ns2 NUd2 n num_runs" << std::endl;
        return 1;
    }
    else{
        // Read in params
        double Ns1 = atof(argv[1]);
        double NUd1 = atof(argv[2]);
        double Ns2 = atof(argv[3]);
        double NUd2 = atof(argv[4]);
        int sample_size = atol(argv[5]);
        int num_samples = atol(argv[6]);

        // Create and seed random number generator
        Random random(std::random_device().operator()());
        
        // Simulate the frequency spectrum
        auto frequency_spectrum = simulate_frequency_spectrum(random, Ns1, NUd1, Ns2, NUd2, sample_size, num_samples);
        // Print the results
        for(int i=1;i<frequency_spectrum.size();++i){
            std::cout << frequency_spectrum[i] << " ";
        }
        std::cout << std::endl;
        return 0;
    }
}

std::vector<double> simulate_frequency_spectrum(Random & random, double Ns1, double NUd1, double Ns2, double NUd2, int sample_size, int num_samples){

    double lambda1;
    if(Ns1 < 1e-10 || NUd1 < 1e-10){ 
        // make sure nothing weird happens when NUd=0
        lambda1 = 0;
    }
    else{
        lambda1 = NUd1/Ns1;
    }
    
    double lambda2;
    if(Ns2 < 1e-10 || NUd2 < 1e-10){ 
        // make sure nothing weird happens when NUd=0
        lambda2 = 0;
    }
    else{
        lambda2 = NUd2/Ns2;
    }
    
    auto fitness_distribution = create_fitness_distribution(Ns1, lambda1, Ns2, lambda2, sample_size);
    auto sample_kx = poisson_distribution(lambda1); // used for drawing the initial sample
    auto sample_ky = poisson_distribution(lambda2);

    // Used for drawing particular lineages to mutate and coalesce
    // Basically a hack so we don't have to construct tons of objects
    std::vector<uniform_int_distribution> draw_index;
    for(int i=0; i<=sample_size; i++){
        draw_index.push_back(uniform_int_distribution(0,i-1));
    }

    // Used to store all the events that are currently possible
    std::vector<Event> event_list;
    event_list.reserve(2*fitness_distribution.size()*fitness_distribution[0].size());

    // vector to store the results
    std::vector<double> frequency_spectrum(sample_size,0);

    // Loop over the desired number of samples
    for(int current_num_samples=0; current_num_samples < num_samples; ++current_num_samples){
        
        // Draw the sample        
        int kx_max = 0; // the maximum populated fitness class
        int ky_max = 0; // the maximum populated fitness class
        
        for(int i=0; i < sample_size; ++i){
            int kx = sample_kx(random);
            if(kx > kx_max){
                extend_fitness_distribution_x(fitness_distribution, kx);
                kx_max = kx;
            }
            
            int ky = sample_ky(random);
            if(ky > ky_max){
                extend_fitness_distribution_y(fitness_distribution, ky);
                ky_max = ky;
            }
            
            fitness_distribution[kx][ky].lineages.push_back(1);
        }
        //std::cout << "Drew the sample!\n";

        while(true) { // Simulate the ancestral history
        
            // Enumerate all the possible events
            event_list.clear();
            double total_rate=0;
            for(int kx=0; kx <= kx_max; ++kx){
                for(int ky=0; ky <= ky_max; ++ky){
                    auto num_lineages = fitness_distribution[kx][ky].lineages.size();
                    if(num_lineages > 0){
                        // Mutation events can only happen in a nonempty class
                        event_list.push_back(Event{total_rate+=fitness_distribution[kx][ky].mutation_rate_x*num_lineages, Event::MUTATION_X, kx, ky});
event_list.push_back(Event{total_rate+=fitness_distribution[kx][ky].mutation_rate_y*num_lineages, Event::MUTATION_Y, kx, ky});
                        //std::cout << "Mutation events for " << num_lineages << " lineages in " << kx << " " << ky << std::endl;
                        //std::cout << fitness_distribution[kx][ky].mutation_rate_x << std::endl;
                        //std::cout << fitness_distribution[kx][ky].mutation_rate_y << std::endl;
   
                        if(num_lineages > 1){
                            // Coalescence events can only happen when there are at least two individuals
                            event_list.push_back(Event{total_rate+=fitness_distribution[kx][ky].coalescence_rate*num_lineages*(num_lineages-1)/2, Event::COALESCENCE, kx, ky}); 
                        } 
                    }                 
                }
            }
            //std::cout << "Enumerated events\n";

            // Randomly choose an event according to its weight
            auto event_ptr = std::lower_bound(event_list.begin(), event_list.end(), total_rate*sample_uniform(random)); 
            // Randomly choose the time of this event
            //double time = sample_exponential(random) / total_rate;
            // Or not randomly if all you want is an average
            double time = 1.0/total_rate;            


            // Update the site-frequency spectrum
            for(int kx=0; kx <= kx_max; ++kx){
                for(int ky=0; ky <= ky_max; ++ky){
                    for(auto i : fitness_distribution[kx][ky].lineages){
                        frequency_spectrum[i] += time;
                    }
                }
            }
            

            // Apply the event
            if(event_ptr->type == Event::MUTATION_X){
                //std::cout << "Mutation x event!\n";
                
                auto kx = event_ptr->kx; 
                auto ky = event_ptr->ky; 
                
                //std::cout << kx << " " << ky << std::endl;
                
                
                // Draw a random lineage from this class
                int i = draw_index[fitness_distribution[kx][ky].lineages.size()](random);
                
                // Send it to the previous class
                fitness_distribution[kx-1][ky].lineages.push_back(fitness_distribution[kx][ky].lineages[i]);
                fitness_distribution[kx][ky].lineages.erase(fitness_distribution[kx][ky].lineages.begin()+i);
                if(fitness_distribution[kx][ky].lineages.empty() && kx==kx_max) 
                    ; //--kx_max; 
            }
            else if(event_ptr->type == Event::MUTATION_Y){
                //std::cout << "Mutation y event!\n";
                
                auto kx = event_ptr->kx; 
                auto ky = event_ptr->ky; 
                //std::cout << kx << " " << ky << std::endl;
                                
                // Draw a random lineage from this class
                int i = draw_index[fitness_distribution[kx][ky].lineages.size()](random);
                
                // Send it to the previous class
                fitness_distribution[kx][ky-1].lineages.push_back(fitness_distribution[kx][ky].lineages[i]);
                fitness_distribution[kx][ky].lineages.erase(fitness_distribution[kx][ky].lineages.begin()+i);
                if(fitness_distribution[kx][ky].lineages.empty() && ky==ky_max) 
                    ; //--ky_max; 
            }
            else if(event_ptr->type == Event::COALESCENCE){
                //std::cout << "Coalescence event!\n";
                
                auto kx = event_ptr->kx; 
                auto ky = event_ptr->ky; 
                
                // Draw a random pair of lineages in this class
                int i = draw_index[fitness_distribution[kx][ky].lineages.size()](random);
                int j = draw_index[fitness_distribution[kx][ky].lineages.size()-1](random);
                if(j>=i) 
                    ++j;
                
                // Coalesce the lineages
                fitness_distribution[kx][ky].lineages[i]+=fitness_distribution[kx][ky].lineages[j];
                if(fitness_distribution[kx][ky].lineages[i] < sample_size){
                    fitness_distribution[kx][ky].lineages.erase(fitness_distribution[kx][ky].lineages.begin()+j);
                }
                else{
                    // All the lineages have coalesced.
                    // Ancestral history is complete!
                    fitness_distribution[kx][ky].lineages.clear();
                    break;
                }
            }
            else{
                    std::cout << "This should never happen!\n";
            }
            //std::cout << "Applied event!\n";
        } 
    }

    // Divide by number of samples to calculate the average
    for(auto & f : frequency_spectrum){
        f/=num_samples;
    }

    return frequency_spectrum;
}

FitnessDistribution create_fitness_distribution(double Ns1, double lambda1, double Ns2, double lambda2, int size_hint){
    FitnessDistribution fitness_distribution(1);
    // zero class
    fitness_distribution[0].push_back(FitnessClass(exp(lambda1+lambda2), 0, 0, size_hint));
    if(lambda2 > 0){
        fitness_distribution[0].push_back(FitnessClass(exp(lambda1+lambda2)/lambda2, 0, Ns2, size_hint));
        extend_fitness_distribution_y(fitness_distribution, (int) (1+4*lambda2));
    }
    if(lambda1 > 0){
        fitness_distribution.push_back(fitness_distribution.back());
        for(auto & fitness_class : fitness_distribution.back()){
            fitness_class.coalescence_rate *= 1.0/lambda1;
            fitness_class.mutation_rate_x += Ns1;
        }
        extend_fitness_distribution_x(fitness_distribution, (int) (1+4*lambda1));
    }
    return fitness_distribution;
}

void extend_fitness_distribution_x(FitnessDistribution & fitness_distribution, int k_max){
    double lambda1 = fitness_distribution[0][0].coalescence_rate/fitness_distribution[1][0].coalescence_rate;
    double Ns1 = fitness_distribution[1][0].mutation_rate_x;
    for(int k = fitness_distribution.size(); k<=k_max; ++k){
        fitness_distribution.push_back(fitness_distribution.back());
        for(auto & fitness_class : fitness_distribution.back()){
            fitness_class.coalescence_rate *= k/lambda1;
            fitness_class.mutation_rate_x += Ns1;
            fitness_class.lineages.clear();
        }
    }
}

void extend_fitness_distribution_y(FitnessDistribution & fitness_distribution, int k_max){
    // infer parameters from fitness distribution 
    double lambda2 = fitness_distribution[0][0].coalescence_rate/fitness_distribution[0][1].coalescence_rate;
    double Ns2 = fitness_distribution[0][1].mutation_rate_y;
    
    // extend in column space 
    for(auto & fitness_distribution_y : fitness_distribution){ // for each row
        // extend columns
        for(int k = fitness_distribution_y.size(); k<=k_max; ++k){       
            fitness_distribution_y.push_back( FitnessClass(fitness_distribution_y.back().coalescence_rate*k/lambda2, fitness_distribution_y.back().mutation_rate_x, Ns2*k, fitness_distribution_y.back().lineages.capacity()) );     
        }
    }
}




