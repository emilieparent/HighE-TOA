"""
A module collecting various routines for calculating pulsation test
test statistics and helper functions.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pulsar/stats.py,v 1.4 2011/07/08 23:05:48 kerrm Exp $

author: M. Kerr <matthew.kerr@gmail.com>
"""

import numpy as np

TWOPI = np.pi*2

def vec(func):
    from numpy import vectorize
    return vectorize(func,doc=func.__doc__)

def to_array(x):
    x = np.asarray(x)
    if len(x.shape)==0: return np.asarray([x])
    return x
    
def from_array(x):
    if (len(x.shape)==1) and (x.shape[0]==1): return x[0]
    return x


def sig2sigma(sig,two_tailed=True):
    """Convert tail probability to "sigma" units, i.e., find the value of the 
       argument for the normal distribution beyond which the integrated tail
       probability is sig.  Note that the default is to interpret this number
       as the two-tailed value, as this is the quantity that goes to 0
       when sig goes to 1.
       
       args
       ----
       sig     the change probability

       kwargs
       ------
       two_tailed [True] interpret sig as two-tailed or one-tailed
                          probability in converting
    """
    from scipy.special import erfc,erfcinv
    from scipy.optimize import fsolve

    sig     = to_array(sig)
    results = np.empty_like(sig)

    if np.any( (sig > 1) | (sig <= 0 ) ):
        raise ValueError('Probability must be between 0 and 1.')

    if not two_tailed: sig *= 2

    def inverfc(x,*args):
        return erfc(x/2**0.5)-args[0]

    for isig,mysig in enumerate(sig):
        if mysig < 1e-120: # approx on asymptotic erfc
            x0 = (-2*np.log(mysig*(np.pi)**0.5))**0.5
            results[isig] = x0 - np.log(x0)/(1+2*x0)
        elif mysig > 1e-15:
            results[isig] = erfcinv(mysig)*2**0.5
        else:
            results[isig] = fsolve(inverfc,[8],(sig,))
    return from_array(results)


def sigma2sig(sigma,two_tailed=True):
    """Convert "sigma" units to chance probability, i.e., return the integral
       of the normal distribution from sigma to infinity, or twice that
       quantity if two_tailed.
       args
       ----
       sigma   argument for the normal distribution

       kwargs
       ------
       two_tailed [True] if True, return twice the integral, else once
    """
    # this appears to handle up to machine precision with no problem
    from scipy.special import erfc
    if two_tailed: return erfc(sigma/2**0.5)
    return 1 - 0.5*erfc(-sigma/2**0.5)

def sigma_trials(sigma,trials):
   # correct a sigmal value for a trials factor
   if sigma < 20:
       p = sigma2sig(sigma)*trials
       if p >= 1: return 0
       return sig2sigma(p)
   else:
      # use an asymptotic expansion -- this needs to be checked!
      return (sigma**2 - 2*np.log(trials))**0.5


def z2m(phases,m=2):
    """ Return the Z^2_m test for each harmonic up to the specified m.
        See de Jager et al. 1989 for definition.    
    """

    phases = np.asarray(phases)*TWOPI #phase in radians
    n = len(phases)

    if n < 5e3:  #faster for 100s to 1000s of phases, but requires ~20x memory of alternative

        s = (np.cos(np.outer(np.arange(1,m+1),phases))).sum(axis=1)**2 +\
            (np.sin(np.outer(np.arange(1,m+1),phases))).sum(axis=1)**2

    else:

        s = (np.asarray([(np.cos(k*phases)).sum() for k in range(1,m+1)]))**2 +\
            (np.asarray([(np.sin(k*phases)).sum() for k in range(1,m+1)]))**2

    return (2./n)*np.cumsum(s)

def z2mw(phases,weights,m=2):
   """ Return the Z^2_m test for each harmonic up to the specified m.

       The user provides a list of weights.  In the case that they are
       well-distributed or assumed to be fixed, the CLT applies and the
       statistic remains calibrated.  Nice!
    """

   phases = np.asarray(phases)*(2*np.pi) #phase in radians

   s = (np.asarray([(np.cos(k*phases)*weights).sum() for k in range(1,m+1)]))**2 +\
       (np.asarray([(np.sin(k*phases)*weights).sum() for k in range(1,m+1)]))**2

   return np.cumsum(s) * (2./(weights**2).sum())

def sf_z2m(ts,m=2):
    """ Return the survival function (chance probability) according to the
        asymptotic calibration for the Z^2_m test.
        
        args
        ----
        ts      result of the Z^2_m test
    """
    from scipy.stats import chi2
    return chi2.sf(ts,2*m)


def em_four(phases,m=2,weights=None):
    """ Return the empirical Fourier coefficients up to the mth harmonic.
        These are derived from the empirical trignometric moments."""
   
    phases = np.asarray(phases)*TWOPI #phase in radians

    n = len(phases) if weights is None else weights.sum()
    weights = 1. if weights is None else weights

    aks = (1./n)*np.asarray([(weights*np.cos(k*phases)).sum() for k in range(1,m+1)])
    bks = (1./n)*np.asarray([(weights*np.sin(k*phases)).sum() for k in range(1,m+1)])

    return aks,bks

def em_lc(coeffs,dom):
    """ Evaluate the light curve at the provided phases (0 to 1) for the
        provided coeffs, e.g., as estimated by em_four."""

    dom = np.asarray(dom)*(2*np.pi)

    aks,bks = coeffs
    rval = np.ones_like(dom)
    for i in range(1,len(aks)+1):
        rval += 2*(aks[i-1]*np.cos(i*dom) + bks[i-1]*np.sin(i*dom))
    return rval

def hm(phases,m=20,c=4):
    """ Calculate the H statistic (de Jager et al. 1989) for given phases.
        H_m = max(Z^2_k - c*(k-1)), 1 <= k <= m
        m == maximum search harmonic
        c == offset for each successive harmonic
    """
    phases = np.asarray(phases)*(2*np.pi) #phase in radians

    s = (np.asarray([(np.cos(k*phases)).sum() for k in range(1,m+1)]))**2 +\
        (np.asarray([(np.sin(k*phases)).sum() for k in range(1,m+1)]))**2

    return ((2./len(phases))*np.cumsum(s) - c*np.arange(0,m)).max()


def hmw(phases,weights,m=20,c=4):
    """ Calculate the H statistic (de Jager et al. 1989) and weight each
        sine/cosine with the weights in the argument.  The distribution
        is corrected such that the CLT still applies, i.e., it maintains
        the same calibration as the unweighted version."""

    phases = np.asarray(phases)*(2*np.pi) #phase in radians

    s = (np.asarray([(weights*np.cos(k*phases)).sum() for k in range(1,m+1)]))**2 +\
        (np.asarray([(weights*np.sin(k*phases)).sum() for k in range(1,m+1)]))**2

    return ( (2./(weights**2).sum()) * np.cumsum(s) - c*np.arange(0,m) ).max()


@vec
def sf_hm(h,m=20,c=4):
    """Return (analytic, asymptotic) survival function (1-F(h))
       for the generalized H-test.
       See M. Kerr dissertation for more details and docstring
       for hm.
    """
    if h < 1e-16: return 1.
    from numpy import exp,arange,log,empty
    from scipy.special import gamma
    fact = lambda x: gamma(x+1)
               
    # first, calculate the integrals of unity for all needed orders
    ints = empty(m)
    for i in range(m):
        sv = i - arange(0,i) # summation vector
        ints[i]  = exp(i*log(h+i*c)-log(fact(i)))
        ints[i] -= (ints[:i]*exp(sv*log(sv*c)-log(fact(sv)))).sum()

    # next, develop the integrals in the power series
    alpha = 0.5*exp(-0.5*c)
    return exp(-0.5*h)*(alpha**arange(0,m)*ints).sum()

@vec
def sf_h20_dj1989(h):
    """Convert the H-test statistic to a chance probability according to the
       formula of de Jager et al. 1989 -- NB the quadratic term is NOT correct."""
    if h <= 23:
        return 0.9999755*np.exp(-0.39802*h)
    if h > 50:
        return 4e-8
    return 1.210597*np.exp(-0.45901*h + 0.0022900*h**2)

def sf_h20_dj2010(h):
    """Use the approximate asymptotic calibration from de Jager et al. 2010."""
    return np.exp(-0.4*h)

def sig2h20(sig):
    """Use approximate (de Jager 2010) relation to invert."""
    return -np.log(sig)/0.4
