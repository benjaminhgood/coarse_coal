from math import exp,log,sqrt,pi
from scipy.integrate import quad
from scipy.optimize import newton
from scipy.special import lambertw

def integrand(x,lam):
    return 2*lam*x*(1-exp(-x))*exp(-x)*exp(-2*lam*exp(-x)+lam*exp(-2*x))

def integrandbylam(x,lam):
    return 2*x*(1-exp(-x))*exp(-x)*exp(-2*lam*exp(-x)+lam*exp(-2*x))

def integrandprime(x,lam):
    return 2*x*(1-exp(-x))*exp(-x)*exp(-2*lam*exp(-x)+lam*exp(-2*x))*(1-lam*exp(-x)*(2-exp(-x)))

def calculate_g(lam):
    return quad(lambda x: lam*integrandbylam(x,lam),0,20)[0]

def forward_integrand(x,lam):
    if x*lam < 1e-08:
        return lam
    else:
        return (1-exp(-lam*x))/x

def calculate_forward_g(lam):
    return quad(lambda x: forward_integrand(x,lam), 0, 1)[0]

def calculate_gbylam(lam):
    return quad(lambda x: integrandbylam(x,lam),0,20)[0]

def calculate_gprime(lam):
    return quad(lambda x: integrandprime(x,lam),0,20+log(lam+1))[0]

def condition(x,var):
    lam = var/x/x
    return x*exp(-lam)-calculate_gprime(lam)-calculate_gbylam(lam)/2

def minimum_Ns(x,arg='mutation'):
    if arg=='variance':
        var = x
        guess = sqrt(2*var/lambertw(2*var)).real
        return newton(lambda x: condition(x,var), guess)
    elif arg=='mutation':
        NUd = x
        guess = NUd/lambertw(NUd).real
        return newton(lambda x: condition(x,NUd*x), guess)

def minimum_t2(variance):
    Nseff = minimum_Ns(variance,arg='variance')
    NUdeff = variance/Nseff
    return calculate_t2(Nseff,NUdeff)    

def variance(t2):
    return newton(lambda x: minimum_t2(x)-t2,1/t2)

def approx_g(lam):
    if lam < 0.5:
        return 3.0/2*lam
    else:
        return 0.5772+log(2*lam)-1.0/4/lam #*sqrt(2*pi*exp(-2))*exp(1.0/log(2*lam)-0.25/lam)
        #t = newton(lambda x: 2*lam*x*(1-x)-1+1.0/log(1.0/x),0.5/lam)
        #return 2*lam*exp(-lam)*sqrt(2*pi/(2*lam*t*(1-2*t)+1.0/(log(1.0/t)**2)))*t*(1-t)*log(1.0/t)*exp(lam*(1-t)**2)

def richard_integrand(x,lam):
    if x < 0.01*max([lam**(-0.5),1.0/lam]):
        return (lam+lam*lam)*x
    else:
        return (1-2*exp(-lam*x)+exp(lam*x*x-2*lam*x))/x

def calculate_richard_g(lam):
    if lam < 1e-10:
        return 0
    else:
        return quad(lambda x: richard_integrand(x,lam),0,1)[0]

def calculate_t2(Ns,NUd):
    lam = NUd/Ns
    return exp(-lam)+(calculate_g(lam))/Ns

def calculate_an(n):
    result = 0
    for i in xrange(1,n):
        result += 2.0/i
    return result

def calculate_new_an(n):
    result = 0
    for i in xrange(1,n):
        subresult = 0
        for k in xrange(i+1,n+1):
            top = 1
            bottom = 1
            for j in xrange(i+1,n+1):
                if j!=k:
                    top *= j*(j-1) 
                    bottom *= (j*(j-1)-k*(k-1))
            print top,bottom 
            subresult += (1.0*top)/bottom
        result += 2*subresult/i    
    return result    
    
def calculate_bn(n):
    return n
    result = 0
    for i in xrange(1,n):
        subresult = 0
        for k in xrange(i+1,n+1):
            top = 1
            bottom = 1
            for j in xrange(i+1,n+1):
                if j!=k:
                    top *= j*(j-1) 
                    bottom *= (j*(j-1)-k*(k-1)) 
            subresult += (1.0*k*(k-1)*top)/bottom
        result += subresult/i    
    return result    

def calculate_ttotn(Ns,NUd,n):
    lam = NUd/Ns
    an = calculate_an(n)
    bn = calculate_bn(n)
    return an*exp(-lam)+bn*calculate_g(lam)/Ns

if __name__=='__main__':
    import numpy
    import pylab
    
    lams = numpy.linspace(0,20,100)
    gs = numpy.array([calculate_g(lam) for lam in lams])
    forward_gs = numpy.array([calculate_forward_g(lam) for lam in lams])
    
    pylab.plot(lams,gs,'b-')
    pylab.plot(lams,forward_gs,'r-')
    pylab.show()
