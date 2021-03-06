"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pulsar/toagen.py,v 1.8 2011/09/22 21:51:32 kerrm Exp $

Calculate TOAs with a variety of methods.

Authors: Paul S. Ray <paul.ray@nrl.navy.mil>
         Matthew Kerr <matthew.kerr@gmail.com>
"""


import sys
import numpy as np
import pylab as pl
import scipy.stats
from scipy.optimize import fmin
from .stats import hm,hmw,sf_hm
from .edf import EDF,find_alignment
from collections import deque


try:
    import fftfit
    import psr_utils
except:
    pass

# other module to deal with weights -- find when the significance stops increasing and use
# that set of photons for fits (via min weight);

# also, might want to use a posterior probability?
# i.e., factoring in the timing template; iterate?  will it converge?
# do not want to do this -- taken care of in likelihood itself

SECSPERDAY = 86400.

class TOAGenerator(object):
    """Manage a data set and set of options to produce the required TOAs from LAT events."""

    def init(self):
        pass

    def __init__(self,data,polyco,template,**kwargs):

        self.data = data; self.polyco = polyco; self.template = template
        self.init()
        self.__dict__.update(**kwargs)
        self.mean_err = 1
        # keep track of errors in phase shift for later analysis
        self.phases = deque()
        self.phase_errs = deque()

    def get_toas(self,binner,use_midpoint=True):
        """ Calculate the TOAs specified by the binner."""

        self.phases.clear(); self.phase_errs.clear()

        # Compute observation duration for each TOA
        toas = [] #np.empty(binner.ntoa)
        err_toas = [] #np.empty(binner.ntoa)
        #tim_strings = ['FORMAT 1']
        tim_strings = []
        self.counter = 0

        for ii,(mjdstart,mjdstop) in enumerate(binner):

            # Compute freq and period at middle of observation
            tmid = (mjdstop + mjdstart)/2.
            pe   = self.polyco.find_entry(tmid)
            freq = self.polyco.eval_spin_freq(tmid)
            period = 1.0/freq

            # Compute phase at start of observation or at midpoint
            phase_time = tmid if use_midpoint else mjdstart
            pe = self.polyco.find_entry(phase_time)
            #polyco_phase0 = pe.eval_phase(phase_time)
            polyco_phase0 = self.polyco.eval_phase(phase_time)
            
            # Select phases
            phases,weights = self.data.toa_data(mjdstart,mjdstop)
            if len(phases) == 0: continue
            
            tau,tau_err,prob = self.get_phase_shift(phases,weights,polyco_phase0)
            self.mean_err = (self.mean_err*ii + tau_err)/(ii+1)
            
            # Prepare a string to write to a .tim file or to send to STDOUT
            toa = phase_time + (tau*period)/SECSPERDAY
            toa_err = tau_err*period*1.0e6
            s = " %s 0.0 %.12f %.2f %s -i LAT -np %d -chanceprob %.2e" % (self.data.obs,toa,toa_err,self.data.tel,len(phases),prob)
            toai = int(np.floor(toa))
            toaf = toa - toai
            L = psr_utils.get_str_princeton_toa(toai,toaf,toa_err,freq=0,dm=0,name=self.data.srcname)
            toas.append('%.12f'%toa)
            err_toas.append('%.2f'%toa_err)
            tim_strings.append(L)
            self.counter += 1

        # Note TOAS in MJD, err_toas in microseconds, tim_strings a line for a FORMAT 1 .tim file
        print("")
        return np.asarray(toas),np.asarray(err_toas),tim_strings

class UnbinnedTOAGenerator(TOAGenerator):

    def init(self):
        self.seeds = np.arange(0.01,0.991,0.01)
        self.good_ephemeris = True
        self.phi0 = self.template.get_location()
        self.prev_peak = self.phi0
        self.phase_drift = False
        self.weights = None
        self.plot_stem = None
        self.display = True

    def __toa_error__(self,val,*args):
        f      = self.__toa_loglikelihood__
        delta  = 0.01
        d2 = (f( [val + delta], *args) - 2*f([val], *args) + f( [val - delta], *args))/(delta)**2
        if d2 < 0:
            print(' *** WARNING! Could not estimate error properly.  Setting to a large value...')
            return 0.2
        return d2**-0.5

    def __toa_loglikelihood__(self,p,*args):
        self.template.set_overall_phase(p[0])
        if args[1] is None:
            return -np.log(self.template(args[0])).sum()
        #return -np.log(1+args[1]*(self.template(args[0],suppress_bg=True)-1)).sum()
        return -np.log(1+args[1]*(self.template(args[0])-1)).sum()

    def get_phase_shift(self,phases,weights,polyco_phase0):

        f   = self.__toa_loglikelihood__
        jump = False
        if (self.good_ephemeris):
            # the ephemeris should be good enough that the TOAs don't drift by more than ~0.1 period
            # this allows a good guess at the TOA to prevent a fit to the wrong peak
            
            seed = [self.prev_peak]
            fit  = fmin(f,seed,args=(phases,weights),disp=0,ftol=1e-9,full_output=True)
            jump = abs(self.prev_peak - fit[0][0])/self.mean_err
            if jump > 10:
                print('  Found a jump, doing a blind search now.')
            jump = jump > 10
            peak_shift = (fit[0][0] - self.phi0)
            if abs(peak_shift) > 0.25:
                self.phase_drift = True
            elif (not jump):
                self.prev_peak = fit[0][0]
                tau     = (peak_shift - polyco_phase0)
                tau_err = self.__toa_error__(fit[0][0],phases,weights)
                if self.display: print('          Peak Shift: %.5f +/- %.5f'%(peak_shift,tau_err))
                self.phases.append(peak_shift)
            if self.plot_stem is not None:
                dom1 = np.linspace(0,1,100)
                cod1 = np.asarray([f([x],phases,weights) for x in dom1])
                dom2 = np.linspace(fit[0][0]-0.04,fit[0][0]+0.04,30)
                cod2 = np.asarray([f([x],phases,weights) for x in dom2])
                pl.figure(10); pl.clf();
                ax1 = pl.gca()
                # calculate coordinates for inset - could be more sophisticated with transAxes
                ax2 = pl.axes([0.1+0.5*(fit[0][0]<0.5),0.15,0.25,0.25])
                for i,(dom,cod,ax) in enumerate(zip([dom1,dom2],[cod1,cod2],[ax1,ax2])):
                    ax.plot(dom,cod)
                    ax.axvline(fit[0][0],color='red')
                    ax.axvline(self.phi0,color='k',ls='-')
                    ax.axvline(fit[0][0]-tau_err,color='red',ls='--')
                    ax.axvline(fit[0][0]+tau_err,color='red',ls='--')
                    if i==1:
                        ax.xaxis.set_major_locator(pl.matplotlib.ticker.MaxNLocator(4))
                        ax.axis([dom[0],dom[-1],cod.min(),cod.min()+5])
                        ax.axhline(cod.min()+0.5,color='blue',ls='--')
                name = ('%3d'%(self.counter+1)).replace(' ','0')
                ax1.set_xlabel('(Relative) Phase')
                ax1.set_ylabel('Negative Log Likelihood')
                pl.savefig('%s%s.png'%(self.plot_stem,name))
            
        if (not self.good_ephemeris) or self.phase_drift or jump:
            # look for the absolute best fit without any prior information

            best_ll,tau = np.inf,0
            seed_vals = [f([x],phases,weights) for x in self.seeds]
            top10     = self.seeds[np.argsort(seed_vals)][:2] #NB change
            for seed in top10:
                fit   = fmin(f,[seed],args=(phases,weights),disp=0,ftol=1e-9,full_output=True)
                if fit[1] < best_ll:
                    best_ll = fit[1]
                    tau = fit[0][0]
                    if jump: self.prev_peak = tau

            tau_err = self.__toa_error__(tau,phases,weights)         
            tau -= (self.phi0 + polyco_phase0)
            if self.display: print('  (Blind) Peak Shift: %.5f +/- %.5f'%(tau+polyco_phase0,tau_err))
            self.phases.append(tau+polyco_phase0)        

        self.phase_errs.append(tau_err)
        return tau,tau_err,sf_hm(hm(phases) if (self.weights is None) else hmw(phases,weights))


class BinnedTOAGenerator(TOAGenerator):

    def init(self):
        self.profile_only = False
        self.nbins        = 32
        self.bins         = np.linspace(0,1,self.nbins + 1)

    def _prof_chisq(self,p):
        """
        Compute the Reduced Chi^2 of a binned profile p.  Each bin should hold
        the number of counts in this bin.
        Returns a tuple containing both the reduce Chi^2 and the probability
        that this value occurred by change.
        """
        # Compute mean rate
        mean = p.mean()
        # Compute chi squared
        chisq = ((p-mean)**2/mean).sum()
        # Return reduced chi squared.  DOF is nbins-1 because we are fitting
        # for the mean rate.
        dof = len(p)-1
        redchi = chisq/dof
        return (redchi,scipy.stats.chisqprob(chisq,dof))

    def _measure_phase(self,profile, template, rotate_prof=True):
        """
        measure_phase(profile, template):
            Call FFTFIT on the profile and template to determine the
                following parameters: shift,eshift,snr,esnr,b,errb,ngood
                (returned as a tuple).  These are defined as in Taylor's
                talk at the Royal Society.
        """
        c,amp,pha = fftfit.cprof(template)
        print("template", template)
        pha1 = pha[0]
        if (rotate_prof):
            pha = np.fmod(pha-np.arange(1,len(pha)+1)*pha1,2.0*np.pi)
        print("profile",len(profile),profile)
        print("amp", amp)
        print("pha ",pha)
        shift,eshift,snr,esnr,b,errb,ngood = fftfit.fftfit(profile,amp,pha)
        if 0:
            c2,amp2,pha2=fftfit.cprof(profile)
            pl.ion()
            pl.figure(1)
            pl.title('Amplitudes')
            pl.plot(amp/amp.sum())
            pl.plot(amp2/amp2.sum())
            pl.figure(2)
            pl.title('Phases')
            pl.plot(pha)
            pl.plot(pha2)
            pl.figure(3)
            #pl.title('Profile')
            pl.plot(profile/profile.sum())
            pl.plot(template/template.sum())
            pl.draw()
            tmp = input()
            pl.close('all')
            print("shift, eshift, snr = ",shift,eshift,snr)
        return shift,eshift,snr,esnr,b,errb,ngood

    def get_phase_shift(self,phases,weights,polyco_phase0):
            
        # Compute phase at start of this block
        if options.profile:
            polyco_phase0 = 0.0

        else:
            phases -= polyco_phase0
            phases[phases < 0.] += 1

        profile = np.histogram(phases,bins=self.bins)[0]
        # This is older version from Matthew Kerr
        #profile = np.histogram(phases,bins=self.bins[:-1])[0]


        if not self.profile_only:
            
            # Compute shift from profile    

            # Try using FFTFIT first
            shift,eshift,snr,esnr,b,errb,ngood = self._measure_phase(profile, self.template, True)
            # tau and tau_err are the predicted phase of the pulse arrival
            tau, tau_err = shift/len(profile), eshift/len(profile)
            # Note: "error" flags are shift = 0.0 and eshift = 999.0

            # If that failed, use a time-domain correlation
            if (np.fabs(shift) < 1e-7 and np.fabs(eshift-999.0) < 1e-7):
                   print(" *** Warning!  Bad return from FFTFIT.  Using PRESTO correlation...", file=sys.stderr)
                   # Not enough structure in the template profile for FFTFIT
                   # so use time-domain correlations instead
                   tau = psr_utils.measure_phase_corr(profile, self.template)
                   # This needs to be changed
                   tau_err = 0.1/len(profile)
            
            #tau = np.random.rand(1)[0];tau_err = np.random.rand(1)[0]*0.01 # testing
            redchi,prob = self._prof_chisq(profile)
            return tau,tau_err,prob


        else:
            # This doesn't work very well since it doesn't pause after the any
            # but the first plot.
            pl.ioff() #does this help the pausing problem?
            pl.figure(1)
            bbins = np.arange(2.0*len(profile))/len(profile)
            pl.bar(bbins,np.concatenate((profile,profile)),width=bbins[1])
            #bins = np.arange(2.0*len(profile)+1)/len(profile)
            #pl.plot(bins,list(profile)+list(profile)+[0.0],linestyle='steps-post')
            pl.grid(1)
            pl.xlabel('Pulse Phase')
            pl.ylabel('Number of Photons')
            #pl.title('Profile')
            pl.show()
            of = open('profile.bestprof',"w")
            print("# Profile generated by polyfold.py", file=of)
            for i,np in zip(list(range(len(profile))),profile):
                print("%d %d" % (i,np), file=of)
            of.close()

            return 0,0,0

class EDFTOAGenerator(TOAGenerator):

    def init(self):
        self.edf = EDF(self.phases)

    def get_phase_shift(self,phases,polyco_phase0):

        peak_shift,sups = find_alignment(self.edf,phases)
        pl.figure(3);pl.clf();pl.plot(sups,marker='.');pl.axvline(50);pl.show()

        tau_err = float(input('Estimate error width in delta_phi:'))
        #tau_err = 0.02
        return peak_shift-polyco_phase0,tau_err,sf_hm(hm(phases))
