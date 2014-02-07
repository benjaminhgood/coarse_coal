import numpy
from math import exp
from scipy.stats import hypergeom

def calculate_ks(n):
    return numpy.arange(1,n)*1.0

def calculate_neutral_sfs(n):
    return 1.0/calculate_ks(n)

def calculate_folded_ks(n):
    ks = calculate_ks(n)
    folded_ks = ks[0:(n+1)/2]
    return folded_ks

def fold_sfs(fs):
    n = len(fs)+1
    folded_fs = (fs + fs[::-1])[0:(n+1)/2]
    if (n-1) % 2 != 0:
        folded_fs[-1] *= 0.5
    return folded_fs/folded_fs[0]

def calculate_folded_theta_multiples(n):
    ks = calculate_ks(n)
    folded_ks = ks[0:(n+1)/2]
    theta_multiples = folded_ks*(n-folded_ks)*1.0/n
    if (n-1) % 2 != 0:
        theta_multiples[-1] = folded_ks[-1]
    return theta_multiples/theta_multiples[0]

def calculate_theta_multiples(n):
    return calculate_ks(n)
    
def calculate_maf_estimate(count,n):
    return min([count,n-count])*1.0/n

def calculate_pi_estimate(count, n):
    return 1.0*count*(n-count)/(n*(n-1)/2)

def calculate_downsampled_estimate(i,n,m):
    downsampled_estimate = numpy.zeros(m+1)
    js = numpy.arange(max([0,m-(n-i)]),min([i,m])+1)
    downsampled_estimate[js] += hypergeom.pmf(js,n,i,m)
    return downsampled_estimate[1:-1]

def calculate_pi(fs):
    n = len(fs)+1
    ks = calculate_ks(n)
    return (ks*(n-ks)*fs).sum()/(n*(n-1)/2)

def calculate_maf(f):
    n = len(f)+1
    ks = numpy.arange(1,n)*1.0
    mafs = numpy.array([min(k,n-k) for k in ks])
    return (mafs*f).sum()/f.sum()/n

def downsample_sfs(fs, m):
    n = len(fs)+1
    downsampled_fs = numpy.zeros(m-1)
    for i in xrange(1,n):
        downsampled_fs += calculate_downsampled_estimate(i,n,m)*fs[i-1]
    return downsampled_fs

def calculate_an(n):
    return (1/calculate_ks(n)).sum()

def calculate_waterson(sfs):
    n = len(sfs)+1
    return sfs.sum()/calculate_an(n)

def calculate_schaeffers_d(sfs):
    n = (len(sfs)+1)*1.0
    a1 = (1.0/numpy.arange(1,n)).sum()
    a2 = numpy.square(1.0/numpy.arange(1,n)).sum()
    b1 = (n+2)/(3*(n-1))
    b2 = 2*(n*n+n+3)/(9*n*(n-1))
    c1 = b1-a1
    c2 = b2-(n+2)/(a1*n)+a2/a1/a1
    e1 = c1/a1
    e2 = c2/(a1*a1+a2)
    pi = calculate_pi(sfs)
    S = sfs.sum()
    d = pi-S/a1
    dmin = -(2/n-1/a1)*S
    return d/dmin
    
def calculate_tajimas_d(sfs):
    n = (len(sfs)+1)*1.0
    a1 = (1.0/numpy.arange(1,n)).sum()
    a2 = numpy.square(1.0/numpy.arange(1,n)).sum()
    b1 = (n+1.0)/(3*(n-1))
    b2 = 2*(n*n+n+3.0)/(9*n*(n-1))
    c1 = b1-1.0/a1
    c2 = b2-(n+2.0)/(a1*n)+a2/a1/a1
    e1 = c1/a1
    e2 = c2/(a1*a1+a2)
    pi = calculate_pi(sfs)
    S = sfs.sum()
    d = pi-S/a1
    norm = S*(e2**0.5)
    return d/norm

def calculate_tajimas_d_sufficient(n,pi,S):
    a1 = (1.0/numpy.arange(1,n)).sum()
    a2 = numpy.square(1.0/numpy.arange(1,n)).sum()
    b1 = (n+1.0)/(3*(n-1))
    b2 = 2*(n*n+n+3.0)/(9*n*(n-1))
    c1 = b1-1.0/a1
    c2 = b2-(n+2.0)/(a1*n)+a2/a1/a1
    e1 = c1/a1
    e2 = c2/(a1*a1+a2)
    d = pi-S/a1
    norm = S*(e2**0.5)
    return d/norm


if __name__=='__main__':
    n = 40
    m = 10
    
    neutral_sfs = calculate_neutral_sfs(n)
    downsampled_sfs = downsample_sfs(neutral_sfs, m)
    print "Original:", neutral_sfs
    print "Downsampled:", downsampled_sfs

