#!/usr/bin/env python
"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pulsar/itemplate.py,v 1.4 2012/01/27 19:10:05 kerrm Exp $

Provide a method for interactively fitting a multi-gaussian template to data.

Authors: Paul S. Ray <paul.ray@nrl.navy.mil>
         Matthew Kerr <matthew.kerr@gmail.com>
"""
import sys
import numpy as np
import pylab as pl
from astropy.io import fits as pyfits
from pint.templates import lctemplate, lcfitters, lcprimitives
from optparse import OptionParser
import subprocess
from pint.scripts.photonphase import main as photonphasemain

def light_curve(phases,weights=None,nbins=25,ec='blue',ls='solid',label=None,axes=None,fignum=1,nmc=100,template=None):
    if axes is None:
        import pylab as pl; pl.figure(fignum); axes = pl.gca()
    bins = np.linspace(0,1,nbins+1)
    bcs = (bins[1:]+bins[:-1])/2
    nph = len(phases)

    cod = axes.hist(phases,bins=bins,weights=weights,density=True,histtype='step',ec=ec)[0]
    
    if weights is None:
        err = (cod*float(nbins)/len(phases))**0.5
    else:
        err = np.empty([nbins,nmc])
        for i in xrange(nmc):
            rweights = weights[np.argsort(np.random.rand(nph))]
            err[:,i] = np.histogram(phases,bins=bins,weights=rweights,normed=True)[0]
        err = np.std(err,axis=1)
    axes.errorbar(bcs,cod,yerr=err,color=ec,capsize=0,ls=' ',marker=' ')

    if (weights is not None):
        bg_level = 1-(weights**2).sum()/weights.sum()
        axes.axhline(bg_level,color=ec)
    else: bg_level = 0

    if template is not None:
        dom = np.linspace(0,1,101)
        axes.plot(dom,template(dom)*(1-bg_level)+bg_level,color='red')
        #axes.plot(dom,template(dom,suppress_bg=(weights is not None))*(1-bg_level)+bg_level,color='red')
    return cod

def write_prof_data(data,outfile='profile.dat'):
    fout = open(outfile,'w')
    for i,d in enumerate(data):
        l = " %3d  %f "%(i,d)
        print(l)
        fout.write(l+'\n')
    print("\nWrote profile data in %s\n"%outfile)
    fout.close()

def has_phase_info(ft1name):
	f 		= pyfits.open(ft1name)
	headcols 	= f['EVENTS'].columns
	has_phase	= False
	if 'PULSE_PHASE' in headcols.names:
		has_phase = True
	else:
		print('Event file does not contain photon phases!')
		print("Do you wish to compute photon phases with PINT's photonphase tool?")
		response = input('Y or N ?  ')
		if 'Y' in response or 'y' in response:
			parfile = input("Provide the name of the parfile for phase prediction:  ")
			if len(parfile) < 1:
				print("Invalid parfile.. \n\nExiting.")
				sys.exit()
			extension = '.'+ft1name.split('.')[-1]
			outfile = ft1name.split(extension)[0] + '_PHASE.fits'
			outplot = ft1name.split(extension)[0] + '_photonfold.png'
			def_args = ['--addphase','--absphase','--barytime','--plot', '--plotfile',outplot,'--outfile',outfile,ft1name,parfile]
			photonphasemain(def_args)
			ft1name = outfile  
			print("Now using file: \n %s\n\n"%outfile)
			has_phase = True
		else:
			print("Not computing photon phases. \n\nExiting.")
			sys.exit()
	f.close()
	return ft1name, has_phase


def get_phases(ft1file,get_weights=False,weightcol='WEIGHT'):
	f = pyfits.open(ft1file)
	ft1file, has_phase = has_phase_info(ft1file)
	if not has_phase:
		print("FT1FILE must contain PULSE_PHASE header entries. \n\nExiting.")
		sys.exit()
	f = pyfits.open(ft1file)
	phases = np.asarray(f['EVENTS'].data.field('PULSE_PHASE'),dtype=float)
	if get_weights:
		weights = np.asarray(f['EVENTS'].data.field(weightcol),dtype=float)
	else: weights = None
	f.close()
	return ft1file,phases,weights


class InteractiveFitter(object):

    def init(self,nbins):
        self.nbins = nbins
        self.weights = None
        self.fignum = 1

    def welcome(self):
        print('Welcome to the interactive unbinned template fitter!')
        print('Displaying the profile... now, we will specify where to put Gaussians.')
        print('For each peak, drag a horizontal line')
        print('         AT THE HIGHEST POINT ON THE PEAK')
        print('         from HALF-MAX to HALF-MAX')
        print('After each drag, the plot will refresh with the current template.')
        print('After all Gaussians are specified, close the plot, and fitting will start.')
        print('(Note -- if using interactively, you will start the fit with do_fit; but do close the plot!')
        
    def __init__(self,inputfn,phases,nbins,Presto,**kwargs):
        self.init(nbins)
        self.fname = inputfn
        self.__dict__.update(**kwargs)
        self.phases = phases
        self.primitives = []
        self.dom = np.linspace(0,1,100)
        self.welcome()
        pl.close(self.fignum)
        self.fig = pl.figure(self.fignum)
        self.ax  = pl.gca()
        self.connect()
        histo = light_curve(self.phases,weights=self.weights,nbins=self.nbins,axes=self.ax)
        self.prof_hist = np.asarray(histo)
        self.outfile = inputfn.replace('.evt','_prof_%dbins.dat'%nbins).replace('.fits','_prof_%dbins.dat'%nbins)
        write_prof_data(self.prof_hist,outfile=self.outfile)
        if Presto:
            pl.close()
        else:
            pl.show()

    def do_fit(self):
        print('Fitting the template with unbinned likelihood...')
        template = lctemplate.LCTemplate(self.primitives)
        fitter   = lcfitters.LCFitter(template,self.phases,weights=self.weights)
        fitter.fit()
        print('Fitting finished!')
        print(fitter)
        print('Overlaying fitted template...')
        self.fig = pl.figure(self.fignum)
        self.ax = pl.gca()
        light_curve(self.phases,weights=self.weights,nbins=self.nbins,axes=self.ax,template=template)
        pl.show()
        self.fitter = fitter

    def connect(self):
        self.cidpress  = self.fig.canvas.mpl_connect('button_press_event',self.on_press)
        self.cidrelese = self.fig.canvas.mpl_connect('button_release_event',self.on_release)

    def on_press(self,event):
        self.x0 = event.xdata
        self.y0 = event.ydata

    def on_release(self,event):
        x1 = event.xdata
        y1 = event.ydata

        fwhm  = x1 - self.x0
        peak  = (y1 + self.y0)/2.
        phase = (x1 + self.x0)/2.

        # just Gaussian for now
        sigma = fwhm/(8 * np.log(2))**0.5
        ampl  = peak * sigma * (2*np.pi)**0.5

        self.primitives += [lcprimitives.LCGaussian(p=[ampl,sigma,phase])]
        template = lctemplate.LCTemplate(self.primitives)
        self.ax.clear()
        light_curve(self.phases,weights=self.weights,nbins=self.nbins,axes=self.ax,template=template)
        pl.draw()

    def write_template(self,outfile):
        if not hasattr(self,'fitter'):
            print('Must do fit first!'); return
        self.fitter.write_template(outfile)


if __name__ == '__main__':

    desc="Read an FT1 file containing PULSE_PHASE info and interactively fit a template."""
    parser=OptionParser(usage=" %prog [options] [FT1_FILENAME]", description=desc)
    parser.add_option('-n','--nbins',type='int',default=50,help="Number of bins to use in phase histogram.")
    parser.add_option('-w','--weights',action='store_true',default=False,help='Use weighted light curve')
    parser.add_option('-c','--weightcol',type='string',default='WEIGHT',help='Column in FT1 file that holds the weight')
    parser.add_option('-p','--prof',type='string',default=None,help='Output name for profile')
    parser.add_option('-m','--min_weight',type='float',default=1e-2,help='Minimum weight to include in fit.')
    
    ## Parse arguments
    (options,args) = parser.parse_args()
    if len(args) < 1:
        parser.print_help()
        raise ValueError('Must provide an FT1 file!')

    inputfn = args[0]
    inputfn, phases,weights = get_phases(inputfn,get_weights=options.weights,weightcol=options.weightcol)

    if options.weights:
        phases = phases[weights > options.min_weight]
        print('%d of %d photons survive weight cut'%(len(phases),len(weights)))
        weights = weights[weights > options.min_weight]

    print('Welcome to the interactive unbinned template fitter!')
    print("What type of template would you like to fit (gauss=Gaussian, kd=Kernel Density, ef [NHARM]=Empirical Fourier, presto=Presto's pygaussfit tool):")
    line = input()
    if line.startswith('gauss'):
        intf = InteractiveFitter(inputfn,phases,options.nbins,Presto=False,weights=weights)
        intf.do_fit()

        if options.prof is not None: out = options.prof
        else:
            out = ''
            out = input('Enter filename for gaussian profile output file, or just hit ENTER to exit...:  ')
        if len(out) > 0:
            print('Writing Gaussian-style template to %s...'%(out))
            intf.write_template(out)

    elif line.startswith('presto'):
        intf = InteractiveFitter(inputfn,phases,options.nbins,Presto=True,weights=weights)
        cmd = "presto_pyGaussfit.py %s"%intf.outfile
        subprocess.call(cmd,shell=True)

    elif line.startswith('kd'):
        dom = np.linspace(0.0,1.0,100)
        prim = lcprimitives.LCKernelDensity(phases=phases)
        template = lctemplate.LCTemplate([prim])
        pl.hist(phases,options.nbins,density=True,histtype='step',edgecolor='k')
        pl.plot(dom,template(dom),color='red')
        pl.title('Kernel Density Template Fit')
        pl.show()
        s = input('Enter a filename here to save template for future use.  Just hit ENTER to skip the step.\n')
        if len(s) > 0:
            prim.to_file(s)
        
    elif line.startswith('ef'):
        dom = np.linspace(0.0,1.0,100)
        if len(line.split()) > 1:
            nharm = int(line.split()[1])
        else:
            nharm = 16
        lcf = lcprimitives.LCEmpiricalFourier(phases=phases,nharm=nharm)
        template = lctemplate.LCTemplate([lcf])
        pl.hist(phases,options.nbins,density=True,histtype='step',edgecolor='k')
        pl.plot(dom,template(dom),color='red')
        pl.title('Empirical Fourier Template with %d harmonics' % (nharm,))
        pl.show()
        s = input('Enter a filename here to save template for future use.  Just hit ENTER to skip the step.\n')
        if len(s) > 0:
            lcf.to_file(s)
    else:
        print("Unrecognized template: ",line)
        sys.exit(1)



    print('Goodbye!')
    
    

