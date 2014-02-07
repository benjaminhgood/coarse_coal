#include <iostream>
#include "evolve.hpp"
#include "phylogeny.hpp"

int single_s_simulation(int argc,char * argv[]);
int uniform_simulation(int argc,char * argv[]);
int finite_sites_simulation(int argc, char * argv[]);
int gamma_simulation(int argc, char * argv[]);
int truncated_exponential_simulation(int argc, char * argv[]);

int main(int argc, char * argv[]){
    return uniform_simulation(argc,argv);
}

int finite_sites_simulation(int argc, char * argv[]){
    if(argc < 9){
        std::cout << "usage: " << argv[0] << " num_samples sample_size N s U L kLi Un" << std::endl;
        return 1;
    }
    else{
        int num_samples = atoi(argv[1]);
        int sample_size = atoi(argv[2]);
        double N = atof(argv[3]);
        double s = atof(argv[4]);
        double U = atof(argv[5]);
        int L = atol(argv[6]);
        double kL0 = atof(argv[7]);
        double Un = atof(argv[8]);
        //Random random = create_random(42);
        FiniteSitesDeltaDFE dfe(U,s,L,kL0);
        Random random = create_random();
        SynonymousTrackedDFE tracked_dfe(N,Un);
        auto results = evolve_population(random,N,dfe,tracked_dfe,num_samples,sample_size);

        auto avg_frequency_spectrum = mean(results.frequency_spectra);
        auto pi = calculate_pi(avg_frequency_spectrum);
        auto Sn = calculate_sn(avg_frequency_spectrum);
        auto sparse_pi = mean(results.pis);
        auto sparse_pi2 = mean_squared(results.pis);

        if(results.avg_fitnesses.size() > 2){
            std::cout << (log(results.avg_fitnesses.back())-log(results.avg_fitnesses[1]))/(results.record_times.back()-results.record_times[1]) << " ";
            std::cout << pi << " " << sparse_pi << " " << sparse_pi2 << " " << Sn << std::endl;
            for(auto fi : avg_frequency_spectrum){
                std::cout << fi << " ";
            }
            std::cout << std::endl;
            for(auto pi_i : results.pis){
                std::cout << pi_i << " "; 
            }
            std::cout << std::endl;
        }
        else{
            std::cout << "Finished!" << std::endl;
        }

        return 0;
    }
}

int single_s_simulation(int argc, char * argv[]){
    if(argc < 7){
        std::cout << "usage: " << argv[0] << " num_samples sample_size pop_size s U Un" << std::endl;
        return 1;
    }
    else{
        int num_samples = atoi(argv[1]);
        int sample_size = atoi(argv[2]);
        double N = atof(argv[3]);
        double s = atof(argv[4]);
        double U = atof(argv[5]);
        double Un = atof(argv[6]);
        //Random random = create_random(43);
        Random random = create_random();
        SynonymousTrackedDFE tracked_dfe(N,Un);
        DeltaDFE dfe(U,s);
        auto results = evolve_population(random,N,dfe,tracked_dfe,num_samples,sample_size);

        auto avg_frequency_spectrum = mean(results.frequency_spectra);
        auto pi = calculate_pi(avg_frequency_spectrum);
        auto Sn = calculate_sn(avg_frequency_spectrum);
        auto sparse_pi = mean(results.pis); 
        auto sparse_pi2 = mean_squared(results.pis);

        if(results.avg_fitnesses.size() > 2){
            std::cout << (log(results.avg_fitnesses.back())-log(results.avg_fitnesses[1]))/(results.record_times.back()-results.record_times[1]) << " ";
            std::cout << pi << " " << sparse_pi << " " << sparse_pi2 << " " << Sn << std::endl;
            for(auto fi : avg_frequency_spectrum){
                std::cout << fi << " ";
            }
            std::cout << std::endl;
            for(auto pi_i : results.pis){
                std::cout << pi_i << " ";
            }
            std::cout << std::endl;
        }
        else{
            std::cout << "Finished!" << std::endl;
        }

        return 0;
    }
}

int uniform_simulation(int argc, char * argv[]){
    if(argc < 8){
        std::cout << "usage: " << argv[0] << " num_samples sample_size pop_size smin smax U Un" << std::endl;
        return 1;
    }
    else{
        int num_samples = atoi(argv[1]);
        int sample_size = atoi(argv[2]);
        double N = atof(argv[3]);
        double smin = atof(argv[4]);
        double smax = atof(argv[5]);
        double U = atof(argv[6]);
        double Un = atof(argv[7]);
        //Random random = create_random(43);
        Random random = create_random();
        SynonymousTrackedDFE tracked_dfe(N,Un);
        UniformDFE dfe(U,smin,smax);
        auto results = evolve_population(random,N,dfe,tracked_dfe,num_samples,sample_size);

        auto avg_frequency_spectrum = mean(results.frequency_spectra);
        auto pi = calculate_pi(avg_frequency_spectrum);
        auto Sn = calculate_sn(avg_frequency_spectrum);
        auto sparse_pi = mean(results.pis);
        auto sparse_pi2 = mean_squared(results.pis);

        if(results.avg_fitnesses.size() > 2){
            std::cout << (log(results.avg_fitnesses.back())-log(results.avg_fitnesses[1]))/(results.record_times.back()-results.record_times[1]) << " ";
            std::cout << pi << " " << sparse_pi << " " << sparse_pi2 << " " << Sn << std::endl;
            for(auto fi : avg_frequency_spectrum){
                std::cout << fi << " ";
            }
            std::cout << std::endl;
        }
        else{
            std::cout << "Finished!" << std::endl;
        }

        return 0;
    }
}


int gamma_simulation(int argc, char * argv[]){
    if(argc < 8){
        std::cout << "usage: " << argv[0] << " num_samples sample_size pop_size savg alpha U Un" << std::endl;
        return 1;
    }
    else{
        int num_samples = atoi(argv[1]);
        int sample_size = atoi(argv[2]);
        double N = atof(argv[3]);
        double savg = atof(argv[4]);
        double alpha = atof(argv[5]);
        double U = atof(argv[6]);
        double Un = atof(argv[7]);
        //Random random = create_random(43);
        Random random = create_random();
        SynonymousTrackedDFE tracked_dfe(N,Un);
        GammaDFE dfe(U,savg,alpha);
        auto results = evolve_population(random,N,dfe,tracked_dfe,num_samples,sample_size);
        auto avg_frequency_spectrum = mean(results.frequency_spectra);
        auto pi = calculate_pi(avg_frequency_spectrum);
        auto Sn = calculate_sn(avg_frequency_spectrum);
        auto sparse_pi = mean(results.pis);
        auto sparse_pi2 = mean_squared(results.pis);

        if(results.avg_fitnesses.size() > 2){
            std::cout << (log(results.avg_fitnesses.back())-log(results.avg_fitnesses[1]))/(results.record_times.back()-results.record_times[1]) << " ";
            std::cout << pi << " " << sparse_pi << " " << sparse_pi2 << " " << Sn << std::endl;
            for(auto fi : avg_frequency_spectrum){
                std::cout << fi << " ";
            }
            std::cout << std::endl;
        }
        else{
            std::cout << "Finished!" << std::endl;
        }

        return 0;
    }
}

int truncated_exponential_simulation(int argc, char * argv[]){
    if(argc < 8){
        std::cout << "usage: " << argv[0] << " num_samples sample_size pop_size smax c U Un" << std::endl;
        return 1;
    }
    else{
        int num_samples = atoi(argv[1]);
        int sample_size = atoi(argv[2]);
        double N = atof(argv[3]);
        double smax = atof(argv[4]);
        double c = atof(argv[5]);
        double U = atof(argv[6]);
        double Un = atof(argv[7]);
        //Random random = create_random(43);
        Random random = create_random();
        SynonymousTrackedDFE tracked_dfe(N,Un);
        TruncatedExponentialDFE dfe(U,smax,c);
        auto results = evolve_population(random,N,dfe,tracked_dfe,num_samples,sample_size);
        auto avg_frequency_spectrum = mean(results.frequency_spectra);
        auto pi = calculate_pi(avg_frequency_spectrum);
        auto Sn = calculate_sn(avg_frequency_spectrum);
        auto sparse_pi = mean(results.pis);
        auto sparse_pi2 = mean_squared(results.pis);

        if(results.avg_fitnesses.size() > 2){
            std::cout << (log(results.avg_fitnesses.back())-log(results.avg_fitnesses[1]))/(results.record_times.back()-results.record_times[1]) << " ";
            std::cout << pi << " " << sparse_pi << " " << sparse_pi2 << " " << Sn << std::endl;
            for(auto fi : avg_frequency_spectrum){
                std::cout << fi << " ";
            }
            std::cout << std::endl;
        }
        else{
            std::cout << "Finished!" << std::endl;
        }

        return 0;
    }
}
