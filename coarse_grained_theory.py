import corrections
import oskar_theory
import coalescent_sim
import sfs_utils
from math import exp

def calculate_variance(Ns,NUd):
    # Calculate the variance in fitness, (Nsigma)^2
    return Ns*Ns*oskar_theory.calculate_var_k(Ns*1.0,NUd*1.0)

def calculate_critical_Ns(variance):
    # This defines the location of the boundary between
    # the interference and background selection regimes
    # in [Ns,(Nsigma)^2] coordinates
    return corrections.minimum_Ns(variance,arg='variance')

def calculate_critical_NUd(Nseff, variance):
    # Assume that at the border,
    # (Nsigma)^2 = NU*Ns
    return variance/Nseff

def calculate_critical_params(variance):
    # This translates the location of the boundary
    # to (Ns, NU) coordinates
    Nseff = calculate_critical_Ns(variance)
    NUdeff = calculate_critical_NUd(Nseff, variance)
    return Nseff, NUdeff

def is_interference_selection(Ns,NUd):
    # True if Ns,NUd falls in the interference selection regime?
    return Ns <= calculate_critical_Ns(calculate_variance(Ns,NUd)) 
    
def calculate_effective_params(Ns,NUd,NUn=0):
    # If (Ns,NUd) falls in the interference selection regime,
    # returns the coarse-grained effective parameters.
    # Otherwise, returns (Ns,NUd,NUn)
    
    Ns = Ns*1.0
    NUd = NUd*1.0
    
    # Calculate the corresponding point on the critical line
    Nseff, NUdeff = calculate_critical_params(calculate_variance(Ns,NUd))

    if Ns >= Nseff:
        # In BGS regime
        return Ns,NUd,NUn
    else:
        # In IS regime
        return Nseff,NUdeff,NUn+NUd-NUdeff

def calculate_t2(Ns,NUd,use_coalescent_sim=False):
    # Calculates mean pairwise coalescent time, T2/N
    # also equal to pi/pi0 from the text
    Nseff,NUdeff,NUneff = calculate_effective_params(Ns,NUd,0)
    if use_coalescent_sim:
        return coalescent_sim.calculate_t2(Nseff,NUdeff,num_samples=100000)
    else:
        return corrections.calculate_t2(Nseff,NUdeff)

def calculate_Q(Ns,NUd,n, num_samples=100000):
    # calculates the average site frequency spectrum 
    # normalized by the number of singletons
    # This is Q_n(i) from the text
    Nseff, NUdeff, NUneff = calculate_effective_params(Ns,NUd,0)
    Q = coalescent_sim.calculate_time_spectrum(Nseff, NUdeff, n, num_samples=num_samples)
    return Q/Q[0] 

def calculate_maf(Ns,NUd,n):
    # Calculates the average minor allele frequency 
    return sfs_utils.calculate_maf(calculate_Q(Ns,NUd,n,num_samples=10000))

def calculate_schaeffers_D(Ns,NUd,n):
    # Calculates the Schaeffers D, 
    # defined as Tajimas D/Dmin    
    return sfs_utils.calculate_schaeffers_D(calculate_Q(Ns,NUd,n,num_samples=10000))
     
if __name__=='__main__':

    NU=50.0
    Ns=1.0
    print "Coarse-graining parameters Ns=%g and NU=%g..." % (Ns,NU)
    maf = calculate_maf(Ns,NU,10)
    print "Success!"