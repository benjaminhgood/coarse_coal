#include "coalescent.hpp"
#include <iostream>

int frequency_spectrum_simulation(int argc, char * argv[]);

int main(int argc, char* argv[]){
    return frequency_spectrum_simulation(argc,argv);
}

int frequency_spectrum_simulation(int argc, char * argv[]){
    if(argc < 6){
        std::cout << "usage: " << argv[0] << " Ns NUd n num_runs" << std::endl;
        return 1;
    }
    else{
        double Ns = atof(argv[1]);
        double NUd = atof(argv[2]);
        double NR = atof(argv[3]);
        int n = atol(argv[4]);
        int num_runs = atol(argv[5]);
        Random random = create_random();
        auto frequency_spectrum = simulate_frequency_spectrum(random,Ns,NUd,NR,n,num_runs);
        for(int i=1;i<frequency_spectrum.size();++i){
            std::cout << frequency_spectrum[i] << " ";
        }
        std::cout << std::endl;
    }
    return 0;
}

