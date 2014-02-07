#include "coalescent.hpp"
#include "lineage.hpp"
#include "event.hpp"
#include <vector>
#include <iostream>
#include <math.h>

typedef std::vector<double> FitnessClassSizes; 

void augment_hks(FitnessClassSizes &, int);

void augment_hks(FitnessClassSizes & hks,int k_max){
    auto lam = hks[1]/hks[0];
    auto k = hks.size()-1;    
    auto hk = hks[k]; 
    hks.reserve(k_max+1);
    while(k < k_max){
        k += 1;
        hk *= lam/k;
        hks.push_back(hk);
    }    
}

class Event{
    // (unfortunately) non-POD data type for an event

    public:
        double location; // used for intrusive random sample
        
        enum Type { MUTATION, COALESCENCE, RECOMBINATION };
        Type type; // whether this event is a coalescence or mutation
        decltype(Sample().begin()) lineage_ptr_ptr, other_lineage_ptr_ptr;
        LineagePtr ancestral_lineage_ptr;
        bool operator<(double p){ return location < p; };
};
typedef std::vector<Event> EventList;

FrequencySpectrum simulate_frequency_spectrum(Random & random, double Ns, double NUd, double NR, int n, int total_runs){

    LineagePool lineage_pool(n*(n-1)/2+n, Lineage());
    
    EventList event_list;
    event_list.reserve(n*(n-1)/2+n+n);

    FrequencySpectrum frequency_spectrum(0.0,n);
    
    Sample sample; 
    sample.reserve(n);

    double lam;
    if(Ns < 1e-10 || NUd < 1e-10){
        lam = 0;
    }
    else{
        lam = NUd/Ns;
    }
    auto sample_k = create_poisson(lam);
   
    auto hks = FitnessClassSizes(2);
    hks[0] = exp(-lam);
    hks[1] = lam*hks[0];
    augment_hks(hks, (int) 4*lam);    

    for(int num_runs=0;num_runs < total_runs;++num_runs){
        
        int k_max = 0;
        // create sample
        for(int i=0;i<n;++i){
            auto k = sample_k(random);
            if(k > k_max) k_max = k;
            auto lineage_ptr = lineage_pool.allocate();
            lineage_ptr->blocks.assign(1, Block{k,(double) k,0,1});
            lineage_ptr->k = k;
            lineage_ptr->num_descendents = 1;
            sample.push_back(lineage_ptr);
        }
        augment_hks(hks,k_max);                
        //std::cout << "Created sample!\n";

        while(sample.size() > 1){
            //std::cout << "About to enumerate events with " << sample.size() << " lineages\n";
            // create list of allowed events;
            event_list.clear(); 
            double total_rate=0;
            for(auto lineage_ptr_ptr = sample.begin(), end=sample.end(); lineage_ptr_ptr != end; ++lineage_ptr_ptr){
                auto k = (*lineage_ptr_ptr)->k;
                
                if(NR > 0){
                    event_list.push_back(Event{total_rate+=NR,Event::RECOMBINATION,lineage_ptr_ptr,lineage_ptr_ptr,*lineage_ptr_ptr});
                }
                //std::cout << "Added recombination event\n";

                if(k > 0){
                    event_list.push_back(Event{total_rate+=Ns*k,Event::MUTATION,lineage_ptr_ptr,lineage_ptr_ptr,*lineage_ptr_ptr});
                }

                //std::cout << "Added mutation event\n";

                for(auto other_lineage_ptr_ptr = lineage_ptr_ptr+1; other_lineage_ptr_ptr != end; ++other_lineage_ptr_ptr){

                    if(k == (*other_lineage_ptr_ptr)->k){
                        auto ancestral_lineage_ptr = lineage_pool.allocate();
                        //std::cout << "Allocated ancestral lineage\n";
                        double rate = calculate_ancestral_lineage(*lineage_ptr_ptr, *other_lineage_ptr_ptr, ancestral_lineage_ptr);
                        //std::cout << "Calculated ancestral lineage\n";
                        //std::cout << rate << "\n";
                        if(rate > 0){
                            event_list.push_back(Event{total_rate+=rate/hks[k], Event::COALESCENCE, lineage_ptr_ptr, other_lineage_ptr_ptr, ancestral_lineage_ptr});
                        }
                    }
                }
            }

            //std::cout << "Enumerated events!\n";


            // Randomly choose an event according to its weight
            auto event_ptr = std::lower_bound(event_list.begin(), event_list.end(), total_rate*sample_uniform(random)); 
            
            // Randomly choose the time of this event
            // double time = sample_exponential(random) / total_rate;
            // Or not randomly if all you want is an average
            double time = 1.0/total_rate;            

            // Update the site-frequency spectrum 
            for(auto & lineage_ptr : sample){
                frequency_spectrum[lineage_ptr->num_descendents] += time;
            }


            // Apply the event
            if(event_ptr->type == Event::MUTATION){

                //std::cout << "Chose mutation event!\n";

                // easy part
                auto lineage_ptr = *(event_ptr->lineage_ptr_ptr);
                
                // now have to find which block it came from
                // and update the k_weights
                auto end = lineage_ptr->blocks.end();
                double weight = sample_uniform(random)*lineage_ptr->k;
                auto block_ptr = std::lower_bound(lineage_ptr->blocks.begin(), lineage_ptr->blocks.end(), Block{0,weight,0,0}, less_cumulative_k);
                /*if(block_ptr == end){
                    for(auto & block : lineage_ptr->blocks){
                        std::cout << block.k << " " << block.cumulative_k << "\n";
                    }
                    std::cout << "Got to end in searching for " << weight << "out of " << lineage_ptr->k << "!\n";
                }*/
                while(block_ptr->k == 0){
                    ++block_ptr;
                    //if(block_ptr == end)
                    //    break;
                }
                /*if(block_ptr == end){
                    for(auto & block : lineage_ptr->blocks){
                        std::cout << block.k << " " << block.cumulative_k << "\n";
                    }
                    std::cout << "Got to end in searching for " << weight << "out of " << lineage_ptr->k << "!\n";
                }*/

                block_ptr->k -= 1;
                for(;block_ptr != end; ++block_ptr){
                    block_ptr->cumulative_k -= 1;
                }
                lineage_ptr->k -= 1;
                
                //int total_k = 0;
                //for(auto & block : lineage_ptr->blocks){
                //    total_k += block.k;
                //}
                //std::cout << "Checked k " << total_k << " " << lineage_ptr->k << std::endl;
                     
            }
            else if(event_ptr->type == Event::RECOMBINATION){
                // the worst one

                //std::cout << "Chose recombination event!\n";

                auto lineage_ptr = *(event_ptr->lineage_ptr_ptr);
                double breakpoint = sample_uniform(random);
                auto breakpoint_block_ptr = std::upper_bound(lineage_ptr->blocks.begin(), lineage_ptr->blocks.end(), Block{0,0,0,breakpoint}, less_breakpoint);
                //std::cout << "Found breakpoint\n";
                int n = breakpoint_block_ptr->k;
                double p = (breakpoint-breakpoint_block_ptr->x_lower)/(breakpoint_block_ptr->x_upper - breakpoint_block_ptr->x_lower);
                int k_below;
                if(n > 0){
                    //std::cout << "About to draw binomial with " << n << " " << p << "\n";
                    k_below = std::binomial_distribution<>(n,p).operator()(random);
                    //std::cout << "Drew binomial\n";
                }
                else{
                    //std::cout << "No mutations here!\n";
                    k_below = 0;
                }
                int k_above = breakpoint_block_ptr->k - k_below;
                if(breakpoint > 0.5){
                    //std::cout << "Breakpoint above\n";
                    //std::cout << k_above << std::endl;        
                    breakpoint_block_ptr->k -= k_above;
                    breakpoint_block_ptr->cumulative_k -= k_above;
                    lineage_ptr->k -= k_above;
                    breakpoint_block_ptr->x_upper = breakpoint;
                    lineage_ptr->blocks.erase(breakpoint_block_ptr+1,lineage_ptr->blocks.end());
                    auto new_k = sample_poisson(random,lam*(1-breakpoint));
                    //std::cout << new_k << std::endl;
                    lineage_ptr->blocks.push_back(Block{new_k,breakpoint_block_ptr->cumulative_k+new_k,breakpoint,1});
                    lineage_ptr->k += new_k;
                }
                else{
                    //std::cout << "Breakpoint below\n";
                    breakpoint_block_ptr->k -= k_below;
                    breakpoint_block_ptr->x_lower = breakpoint;
                    auto new_lineage_ptr = lineage_pool.allocate();
                    new_lineage_ptr->num_descendents = lineage_ptr->num_descendents;
                    new_lineage_ptr->blocks.assign(1,Block{sample_poisson(random,lam*breakpoint),0,0,breakpoint});
                    new_lineage_ptr->blocks.insert(new_lineage_ptr->blocks.end(), breakpoint_block_ptr, lineage_ptr->blocks.end());
                    event_ptr->lineage_ptr_ptr->swap(new_lineage_ptr);
                }

                int cumulative_k = 0;
                for(auto & block : (*event_ptr->lineage_ptr_ptr)->blocks){
                    cumulative_k += block.k;    
                    block.cumulative_k = cumulative_k;
                }
                (*event_ptr->lineage_ptr_ptr)->k = cumulative_k;
          
                //int total_k = 0;
                //for(auto & block : (*event_ptr->lineage_ptr_ptr)->blocks){
                //    total_k += block.k;
                //}
                //std::cout << "Checked k " << total_k << " " << (*event_ptr->lineage_ptr_ptr)->k << std::endl;
                

            }
            else if(event_ptr->type == Event::COALESCENCE){

            //std::cout << "Chose coalescent event!\n";

                // pretty easy
                event_ptr->lineage_ptr_ptr->swap(event_ptr->ancestral_lineage_ptr);
                sample.erase(event_ptr->other_lineage_ptr_ptr);
            }
            else{

                std::cout << "Should never be here!\n";
            }
            //std::cout << "Applied event!\n";
        }
        sample.clear();
    }
    frequency_spectrum /= total_runs; 
    return frequency_spectrum;
}
