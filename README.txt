LICENSE: The code located here is provided as-is, with no warranty, etc. (It's under GPL v2.) 

REQUIREMENTS: Python (v2.7) with the NumPy (v1.7) and SciPy (v0.12) libraries, and a C++ compiler that supports the C++11 standard (e.g., GNU's g++). 

INSTALLATION: (1) Compile the "structured_coalescent.cpp" source file and name the program "simulate_time_spectrum". (2) Run "python coarse_grained_theory.py" to ensure that it outputs "Success!". 

USAGE: (1) Structured coalescent simulations can be called from the command line:

./simulate_time_spectrum Ns NU n num_runs

where n is the sample size and num_runs is the number of genealogies to average over. This program outputs the average site frequency spectrum P_n(i) in units of the scaled neutral mutation rate. For convenience, we have also implemented a Python interface to this command line program in "coalescent_sim.py".  

(2) Coarse-grained coalescent predictions are accessible via our Python interface. We have provided an "example.py" script that demonstrates how to call some of these methods. Predictions for asexual populations are implemented in "coarse_grained_theory.py". The coarse-grained parameters can be obtained from the "calculate_effective_params" function. One can also obtain the coarse-grained predictions directly, without this intermediate step, using the methods "calculate_t2", "calculate_Q", "calculate_maf", or "calculate_schaeffers_D". The linkage block approximation is implemented in "recombination_theory.py". The effective asexual parameters can be obtained from the "calculate_effective_params" function. For convenience, we have also defined recombinant versions of "calculate_t2", "calculate_Q", "calculate_maf", and "calculate_tajimas_D". 

AUXILLIARY PROGRAMS: For completeness, we have also included C++ sources for the auxilliary programs utilized in the main text: (1) asexual forward time simulations, (2) recombining forward-time simulations, (3) the Bolthausen-Sznitman coalescent, and (4) the recombining structured coalescent. 

JOURNAL REFERENCE: Good BH, Walczak AM, Neher RA, and Desai MM, "Genetic diversity in the interference selection limit," PLoS Genetics (2013). 
 