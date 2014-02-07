#ifndef LINEAGE_HPP
#define LINEAGE_HPP

#include<vector>
#include "object_pool.hpp"

class Block{
    public:
    //POD data type to keep track of the number of mutations in each segment
    int k;
    double cumulative_k;
    double x_lower;
    double x_upper;
};

auto less_breakpoint = [](Block const & block1, Block const & block2){ return block1.x_upper < block2.x_upper; };
auto less_cumulative_k = [](Block const & block1, Block const & block2){ /*std::cout << block1.cumulative_k << " " << block2.cumulative_k << std::endl;*/ return (block1.cumulative_k < block2.cumulative_k); };

class Lineage{
    public:
        std::vector<Block> blocks; // should be roughly T2*R of them
        int k;
        int num_descendents;        
};

typedef SharedObjectPool<Lineage> LineagePool;
typedef decltype(LineagePool().allocate()) LineagePtr; 
typedef std::vector<LineagePtr> Sample;

inline double calculate_ancestral_lineage(LineagePtr const & lineage1_ptr, LineagePtr const & lineage2_ptr, LineagePtr & ancestral_lineage_ptr){

    // by assumption, once the lineages have gotten to this point, 
    // we have already checked that k1 = k2

    ancestral_lineage_ptr->blocks.clear();
    ancestral_lineage_ptr->k = lineage1_ptr->k;
    ancestral_lineage_ptr->num_descendents = lineage1_ptr->num_descendents+lineage2_ptr->num_descendents;
 
    auto block1_ptr = lineage1_ptr->blocks.begin();
    auto block2_ptr = lineage2_ptr->blocks.begin();

    auto end1 = lineage1_ptr->blocks.end();     
    auto end2 = lineage2_ptr->blocks.end();     

    auto block1 = *block1_ptr; 
    auto block2 = *block2_ptr; 

    double rate=exp(gammaln(block1.k+1))/pow(block1.x_upper-block1.x_lower, block1.k)*exp(gammaln(block2.k+1))/pow(block2.x_upper-block2.x_lower, block2.k);
    //std::cout << rate << std::endl;
    int cumulative_k = 0;
    Block block{0,0,0,0};
    
    /*for(auto & block : lineage1_ptr->blocks){
        std::cout << "(" << block.x_lower << ", " << block.x_upper << "): " << block.k << std::endl;  
    }

    for(auto & block : lineage2_ptr->blocks){
        std::cout << "(" << block.x_lower << ", " << block.x_upper << "): " << block.k << std::endl;  
    }*/

    while(block1_ptr != end1 || block2_ptr != end2){

        if(block1.x_upper > block2.x_upper){
            block.x_upper = block2.x_upper; 
            block.k = block2.k;
            block1.k -= block2.k;
            ++block2_ptr;
            if(block2_ptr != end2){
                block2 = *block2_ptr;
                rate*=exp(gammaln(block2.k+1))/pow(block2.x_upper-block2.x_lower, block2.k);
            } 
        }
        else if(block1.x_upper == block2.x_upper){
            block.x_upper = block2.x_upper;
            if(block1.k == block2.k){
                block.k = block2.k;
            }
            else{
                return -1;
            }

            ++block1_ptr;
            if(block1_ptr != end1){
                block1 = *block1_ptr;
                rate*=exp(gammaln(block1.k+1))/pow(block1.x_upper-block1.x_lower, block1.k);
                //std::cout << rate << std::endl;
            }
            
            ++block2_ptr;
            if(block2_ptr != end2){
                block2 = *block2_ptr;
                rate*=exp(gammaln(block2.k+1))/pow(block2.x_upper-block2.x_lower, block2.k);
                //std::cout << rate << std::endl;
            }     
        }
        else{
            block.x_upper = block1.x_upper;
            block.k = block1.k;
            block2.k -= block1.k;
            ++block1_ptr;
            if(block1_ptr != end1){
                block1 = *block1_ptr;
                rate*=exp(gammaln(block1.k+1))/pow(block1.x_upper-block1.x_lower,block1.k);
                //std::cout << rate;
            }  
        }
        if(block.k < 0){
            return -1;
        }
        cumulative_k += block.k;
        block.cumulative_k = cumulative_k;
        rate *= pow(block.x_upper-block.x_lower, block.k)/exp(gammaln(block.k+1));
        //std::cout << rate << std::endl;
        ancestral_lineage_ptr->blocks.push_back(block);
        block.x_lower = block.x_upper;
    }    
    ancestral_lineage_ptr->k = cumulative_k;
    rate *= exp(-gammaln(ancestral_lineage_ptr->k+1.0));
    /*for(auto & block : ancestral_lineage_ptr->blocks){
        std::cout << "(" << block.x_lower << ", " << block.x_upper << "): " << block.k << std::endl;  
    }*/
    //std::cout << " " << ancestral_lineage_ptr->k << " " << rate << std::endl;
    return rate;
}


#endif
