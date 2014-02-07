import coarse_grained_theory
import recombination_theory

# Asexual population

Ns = 10
NU = 300
NR = 0
n = 10

# Calculate coarse-grained parameters
Nseff,NUeff,NUneff = coarse_grained_theory.calculate_effective_params(Ns,NU)

print "Original parameters:", Ns,NU,NR
print "Coarse-grained parameters:", Nseff,NUeff,0


# Calculate the shape of the SFS directly
Q = coarse_grained_theory.calculate_Q(Ns,NU,n)

print "Qn(i):", Q

# Recombining population

Ns = 10
NU = 300
NR = 300

# Calculate coarse-grained parameters
Nseff,NUeff = recombination_theory.calculate_effective_params(Ns,NU,NR)

print "Original parameters:", Ns,NU,NR
print "Coarse-grained parameters:", Nseff,NUeff,0

# Calculate the shape of the SFS directly
Q = recombination_theory.calculate_Q(Ns,NU,NR,n)

print "Qn(i):", Q