import subprocess
import numpy
from math import log,log10,exp
import sys

disabled=False

if sys.platform.startswith("linux") or sys.platform.startswith("darwin"):
    LINUX = True
else:
    LINUX = False

def calculate_time_spectrum(Ns,NUd,n,num_samples=1000):
    if disabled:
        return numpy.zeros(n-1)
    if LINUX:
        exec_name = "./simulate_time_spectrum"
    else:
        exec_name = "./simulate_time_spectrum.exe"
    sys.stderr.write("Calculating time spectrum for parameters Ns=%g, NU=%g, n=%d\n" % (Ns, NUd, n))
    args = [exec_name,str(Ns),str(NUd),str(n),str(num_samples)]
    output = subprocess.Popen(args,shell=False,stdout=subprocess.PIPE).communicate()[0]
    fs = numpy.array([float(item) for item in output.split()])
    return fs

    return calculate_recombination_time_spectrum(Ns,NUd,0,n,num_samples)
        
def calculate_ttotn(Ns,NUd,n,num_samples=1000):
    return calculate_time_spectrum(Ns,NUd,n,num_samples).sum()

def calculate_t2(Ns,NUd,num_samples=1000):
    return calculate_ttotn(Ns,NUd,2,num_samples)/2