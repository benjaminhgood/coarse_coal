#include <iostream>
#include <vector>
#include <valarray>
#include <cmath>
#include <random>
#include <algorithm>

using namespace std;
typedef std::mt19937 Random;
std::uniform_real_distribution<> sample_uniform(0,1);
std::exponential_distribution<> sample_exponential(1);

inline Random create_random(unsigned int seed=0, bool print_seed=false){
    if(seed==0){
        std::random_device rd;
        seed = rd();
    }
    if(print_seed){
        std::cout << seed << std::endl;
    }
    return Random(seed);
}


inline int draw_num_merged(Random & random, int b){
    return ceil(1.0/(1-(1-1.0/b)*sample_uniform(random)));
}

inline std::valarray<double> simulate_time_spectrum(Random & random, int n, int total_runs){

    std::valarray<double> time_spectrum(0.0,n);
    // zero component is num samples, rest is time spectrum

    std::vector<int> sample(n), new_sample(n);
    for(int num_runs=0;num_runs < total_runs;++num_runs){
        // reset the sample
        sample.assign(n,1);
        new_sample.clear();

        // draw the coalescent tree
        int b = sample.size();
        while(b > 1){
            double time = sample_exponential(random)/(b-1);
            int new_weight = 0;
            for(int i=0,j=0,k=draw_num_merged(random,b);i<b;++i){
                int weight = sample[i];
                time_spectrum[weight]+=time;
           
                if(j < k && (b-i)*sample_uniform(random) < k-j){
                    new_weight += weight;
                    ++j;
                }
                else{
                    new_sample.push_back(weight);
                }                
            }
            new_sample.push_back(new_weight);
            sample.swap(new_sample);
            new_sample.clear();
            b = sample.size();
        }
        time_spectrum[0]+=1.0;
    } 
    return time_spectrum;
}

int main(int argc, char* argv[]){
    if(argc < 3){
        std::cout << "usage: " << argv[0] << " n num_runs" << std::endl;
        return 1;
    }
    else{
        int n = atol(argv[1]);
        int num_runs = atol(argv[2]);
        Random random = create_random();
        /*valarray<double> time_spectrum(0.0,n);
        for(int i=0;i<num_runs;++i){
            time_spectrum[draw_num_merged(random,n)-1]+=1.0;
        }*/
        auto time_spectrum = simulate_time_spectrum(random,n,num_runs);
        time_spectrum /= num_runs; //time_spectrum[0];
        for(int i=1;i<time_spectrum.size();++i){
            cout << time_spectrum[i] << ", ";
        }
        cout << endl;
    }
    return 0;
}

