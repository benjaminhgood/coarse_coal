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
        double mutation_rate; // rate at which lineages leave this class due to a mutation event
        std::vector<int> lineages; // number of present day descendents of each lineage

        FitnessClass(double coalescence_rate=0, double mutation_rate=0, int size_hint = 1): coalescence_rate(coalescence_rate), mutation_rate(mutation_rate) { lineages.reserve(size_hint); }
};

typedef std::vector<FitnessClass> FitnessDistribution;

class Event{
    public:
        enum EventType { MUTATION, COALESCENCE};
        double weight; // used for intrusive search
        EventType type;
        decltype(FitnessDistribution().begin()) k_ptr; // the fitness class where this event occurs
        bool operator<(const double p) const { return weight < p; } // used for intrusive search
};

// The primary function in this file. Simulates the synonymous site frequency spectrum.
std::vector<double> simulate_frequency_spectrum(Random & random, double Ns, double NUd, int sample_size, int num_samples);
// ...
double simulate_ttotn_ratio(Random & random, double Ns, double NUd, int sample_size, int num_samples);

// Initializes the population fitness distribution according to the supplied parameters
// Doesn't create all the fitness classes initially, but only up to 4*lambda+1
FitnessDistribution create_fitness_distribution(double Ns, double lambda, int size_hint);

// Extends the fitness distribution to the new k_max if necessary
void extend_fitness_distribution(FitnessDistribution & fitness_distribution, int k_max);

// Main function. Reads in parameters, calls simulate_frequency_spectrum, and prints results.
int main(int argc, char* argv[]){
    if(argc < 5){
        std::cout << "usage: " << argv[0] << " Ns NUd n num_runs" << std::endl;
        return 1;
    }
    else{
        // Read in params
        double Ns = atof(argv[1]);
        double NUd = atof(argv[2]);
        int sample_size = atol(argv[3]);
        int num_samples = atol(argv[4]);

        // Create and seed random number generator
        Random random(std::random_device().operator()());
      
        
        // Simulate the frequency spectrum
        auto frequency_spectrum = simulate_frequency_spectrum(random, Ns, NUd, sample_size, num_samples);
        // Print the results
        for(int i=1;i<frequency_spectrum.size();++i){
            std::cout << frequency_spectrum[i] << " ";
        }
        std::cout << std::endl;
    
        // Simulate the variance in ttotn
        //std::cout << simulate_ttotn_ratio(random, Ns, NUd, sample_size, num_samples) << std::endl;
 
    }
    return 0;
}

std::vector<double> simulate_frequency_spectrum(Random & random, double Ns, double NUd, int sample_size, int num_samples){

    double lambda;
    if(Ns < 1e-10 || NUd < 1e-10){ 
        // make sure nothing weird happens when NUd=0
        lambda = 0;
    }
    else{
        lambda = NUd/Ns;
    }
    auto fitness_distribution = create_fitness_distribution(Ns, lambda, sample_size);
    auto sample_k = poisson_distribution(lambda); // used for drawing the initial sample

    // Used for drawing particular lineages to mutate and coalesce
    // Basically a hack so we don't have to construct tons of objects
    std::vector<uniform_int_distribution> draw_index;
    for(int i=0; i<=sample_size; i++){
        draw_index.push_back(uniform_int_distribution(0,i-1));
    }

    // Used to store all the events that are currently possible
    std::vector<Event> event_list;
    event_list.reserve(2*fitness_distribution.size());

    // vector to store the results
    std::vector<double> frequency_spectrum(sample_size,0);

    // Loop over the desired number of samples
    for(int current_num_samples=0; current_num_samples < num_samples; ++current_num_samples){
        
        // Draw the sample        
        int k_max = 0; // the maximum populated fitness class
        for(int i=0; i < sample_size; ++i){
            int k = sample_k(random);
            if(k > k_max){
                extend_fitness_distribution(fitness_distribution, k);
                k_max = k;
            }
            fitness_distribution[k].lineages.push_back(1);
        }

        auto kmin_ptr = fitness_distribution.begin(); // Minimum populated fitness class
        auto kmax_ptr = kmin_ptr + k_max; // Maximum populated fitness class

        while(true) { // Simulate the ancestral history
        
            // Enumerate all the possible events
            event_list.clear();
            double total_rate=0;
            for(auto k_ptr = kmin_ptr; k_ptr <= kmax_ptr; ++k_ptr){           
                auto num_lineages = k_ptr->lineages.size();
                if(num_lineages > 0){
                    // Mutation events can only happen in a nonempty class
                    event_list.push_back(Event{total_rate+=k_ptr->mutation_rate*num_lineages, Event::MUTATION, k_ptr});
                    if(num_lineages > 1){
                        // Coalescence events can only happen when there are at least two individuals
                        event_list.push_back(Event{total_rate+=k_ptr->coalescence_rate*num_lineages*(num_lineages-1)/2, Event::COALESCENCE, k_ptr}); 
                    } 
                }                 
            }

            // Randomly choose an event according to its weight
            auto event_ptr = std::lower_bound(event_list.begin(), event_list.end(), total_rate*sample_uniform(random)); 
            // Randomly choose the time of this event
            //double time = sample_exponential(random) / total_rate;
            // Or not randomly if all you want is an average
            double time = 1.0/total_rate;            


            // Update the site-frequency spectrum 
            for(auto k_ptr = kmin_ptr; k_ptr <= kmax_ptr; ++k_ptr){
                for(auto i : k_ptr->lineages){
                    frequency_spectrum[i] += time;
                }
            }

            // Apply the event
            if(event_ptr->type == Event::MUTATION){
                auto k_ptr = event_ptr->k_ptr;
                
                // Draw a random lineage from this class
                int i = draw_index[k_ptr->lineages.size()](random);
                
                // Send it to the previous class
                (k_ptr-1)->lineages.push_back(k_ptr->lineages[i]);
                k_ptr->lineages.erase(k_ptr->lineages.begin()+i);
                if(k_ptr->lineages.empty() && k_ptr==kmax_ptr) 
                    --kmax_ptr; 
            }
            else if(event_ptr->type == Event::COALESCENCE){
                auto k_ptr = event_ptr->k_ptr;
                
                // Draw a random pair of lineages in this class
                int i = draw_index[k_ptr->lineages.size()](random);
                int j = draw_index[k_ptr->lineages.size()-1](random);
                if(j>=i) 
                    ++j;
                
                // Coalesce the lineages
                k_ptr->lineages[i]+=k_ptr->lineages[j];
                if(k_ptr->lineages[i] < sample_size){
                    k_ptr->lineages.erase(k_ptr->lineages.begin()+j);
                }
                else{
                    // All the lineages have coalesced.
                    // Ancestral history is complete!
                    k_ptr->lineages.clear();
                    break;
                }
            }
            else{
                    std::cout << "This should never happen!\n";
            }
        } 
    }

    // Divide by number of samples to calculate the average
    for(auto & f : frequency_spectrum){
        f/=num_samples;
    }

    return frequency_spectrum;
}

double simulate_ttotn_ratio(Random & random, double Ns, double NUd, int sample_size, int num_samples){

    double lambda;
    if(Ns < 1e-10 || NUd < 1e-10){ 
        // make sure nothing weird happens when NUd=0
        lambda = 0;
    }
    else{
        lambda = NUd/Ns;
    }
    auto fitness_distribution = create_fitness_distribution(Ns, lambda, sample_size);
    auto sample_k = poisson_distribution(lambda); // used for drawing the initial sample

    // Used for drawing particular lineages to mutate and coalesce
    // Basically a hack so we don't have to construct tons of objects
    std::vector<uniform_int_distribution> draw_index;
    for(int i=0; i<=sample_size; i++){
        draw_index.push_back(uniform_int_distribution(0,i-1));
    }

    // Used to store all the events that are currently possible
    std::vector<Event> event_list;
    event_list.reserve(2*fitness_distribution.size());

    // vector to store the results
    std::vector<double> frequency_spectrum(sample_size,0);
    double ttotn=0;
    double ttotn2=0;

    // Loop over the desired number of samples
    for(int current_num_samples=0; current_num_samples < num_samples; ++current_num_samples){
        
        // Draw the sample        
        int k_max = 0; // the maximum populated fitness class
        for(int i=0; i < sample_size; ++i){
            int k = sample_k(random);
            if(k > k_max){
                extend_fitness_distribution(fitness_distribution, k);
                k_max = k;
            }
            fitness_distribution[k].lineages.push_back(1);
        }

        auto kmin_ptr = fitness_distribution.begin(); // Minimum populated fitness class
        auto kmax_ptr = kmin_ptr + k_max; // Maximum populated fitness class

        while(true) { // Simulate the ancestral history
        
            // Enumerate all the possible events
            event_list.clear();
            double total_rate=0;
            for(auto k_ptr = kmin_ptr; k_ptr <= kmax_ptr; ++k_ptr){           
                auto num_lineages = k_ptr->lineages.size();
                if(num_lineages > 0){
                    // Mutation events can only happen in a nonempty class
                    event_list.push_back(Event{total_rate+=k_ptr->mutation_rate*num_lineages, Event::MUTATION, k_ptr});
                    if(num_lineages > 1){
                        // Coalescence events can only happen when there are at least two individuals
                        event_list.push_back(Event{total_rate+=k_ptr->coalescence_rate*num_lineages*(num_lineages-1)/2, Event::COALESCENCE, k_ptr}); 
                    } 
                }                 
            }

            // Randomly choose an event according to its weight
            auto event_ptr = std::lower_bound(event_list.begin(), event_list.end(), total_rate*sample_uniform(random)); 
            // Randomly choose the time of this event
            double time = sample_exponential(random) / total_rate;
            

            // Update the site-frequency spectrum 
            for(auto k_ptr = kmin_ptr; k_ptr <= kmax_ptr; ++k_ptr){
                for(auto i : k_ptr->lineages){
                    frequency_spectrum[i] += time;
                }
            }

            // Apply the event
            if(event_ptr->type == Event::MUTATION){
                auto k_ptr = event_ptr->k_ptr;
                
                // Draw a random lineage from this class
                int i = draw_index[k_ptr->lineages.size()](random);
                
                // Send it to the previous class
                (k_ptr-1)->lineages.push_back(k_ptr->lineages[i]);
                k_ptr->lineages.erase(k_ptr->lineages.begin()+i);
                if(k_ptr->lineages.empty() && k_ptr==kmax_ptr) 
                    --kmax_ptr; 
            }
            else if(event_ptr->type == Event::COALESCENCE){
                auto k_ptr = event_ptr->k_ptr;
                
                // Draw a random pair of lineages in this class
                int i = draw_index[k_ptr->lineages.size()](random);
                int j = draw_index[k_ptr->lineages.size()-1](random);
                if(j>=i) 
                    ++j;
                
                // Coalesce the lineages
                k_ptr->lineages[i]+=k_ptr->lineages[j];
                if(k_ptr->lineages[i] < sample_size){
                    k_ptr->lineages.erase(k_ptr->lineages.begin()+j);
                }
                else{
                    // All the lineages have coalesced.
                    // Ancestral history is complete!
                    k_ptr->lineages.clear();
                    break;
                }
            }
            else{
                    std::cout << "This should never happen!\n";
            }
        }

        double current_ttotn = 0;
        for(auto & f : frequency_spectrum){
            current_ttotn += f;
            f=0;
        }
        ttotn += current_ttotn;
        ttotn2 += current_ttotn*current_ttotn;
    }

    // Divide by number of samples to calculate the average
    ttotn = ttotn / num_samples;
    ttotn2 = ttotn2 / num_samples;
    return ttotn2/(ttotn*ttotn)-1;
    
}

FitnessDistribution create_fitness_distribution(double Ns, double lambda, int size_hint){
    FitnessDistribution fitness_distribution;
    fitness_distribution.push_back(FitnessClass(exp(lambda), 0, size_hint));
    if(lambda > 0){
        fitness_distribution.push_back(FitnessClass(exp(lambda)/lambda, Ns, size_hint));
        extend_fitness_distribution(fitness_distribution, (int) (1+4*lambda));
    }
    return fitness_distribution;
}

void extend_fitness_distribution(FitnessDistribution & fitness_distribution, int k_max){
    double lambda = fitness_distribution[0].coalescence_rate/fitness_distribution[1].coalescence_rate;
    double Ns = fitness_distribution[1].mutation_rate;
    for(int k = fitness_distribution.size(); k<=k_max; ++k){
        fitness_distribution.push_back( FitnessClass(fitness_distribution.back().coalescence_rate*k/lambda, Ns*k, fitness_distribution.back().lineages.capacity()) );     
    }
}




