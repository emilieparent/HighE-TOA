#!/usr/bin/env python
# Program: upolyfold.py
# Authors: Paul S. Ray <paul.ray@nrl.navy.mil>
#          Matthew Kerr <matthew.kerr@gmail.com>
# Description:
# Reads an FT1 file with geocentered event times and
# folds according to a polyco.dat file generated by tempo2 at the geocenter site ("coe")
# Then, fits TOAs and outputs them in tempo2 format.
#
# Modifications, 2022-02: Emilie Parent <parent@ice.csic.es>
# Can now handle barycentered data from X-ray instruments
# Also option to calculate photon phase and save into new fits file (PINT's photonphase)

from optparse import OptionParser
import sys
import numpy as np
from astropy.io import fits as pyfits

from uw.pulsar.phasedata import PhaseData
from uw.pulsar.toabinner import TOABinner,UniformLogBinner
from uw.pulsar.toagen import UnbinnedTOAGenerator
from pint.templates import lctemplate
from pint import fits_utils, models, polycos
from pint.scripts.photonphase import main as photonphasemain
from pint.observatory import get_observatory

def start_message(ft1name):
	print(" ")
	print("       +++++++++++++++++++    ")
	print("   Maximum Likelihood TOA extractor ")
	print("       +++++++++++++++++++    \n")
	print(" Extracting TOAs from file:    ")
	print(" %s"%ft1name)

def has_phase_info(ft1file):
	f 		= pyfits.open(ft1file)
	headcols 	= f['EVENTS'].columns
	has_phase	= False
	if 'PULSE_PHASE' in headcols.names:
		has_phase = True
	else:
		print('Warning: event file does not contain photon phases.')
	f.close()
	return has_phase

def get_time_sys(ft1file):
	f 	= pyfits.open(ft1file)
	ft1hdr  = f['EVENTS'].header
	timesys = ft1hdr['TIMESYS']
	tel 	= ft1hdr['TELESCOP'].lower()
	if timesys=='TDB':
		site = 'barycenter'
	else:
		try:
			site = get_observatory(tel).name
		except:
			site = 'geocenter'

	f.close()
	return timesys,site
	
	


desc="""Read an FT1 file containing GEOCENTERED events (typically generated by gtbary) and fold according to a polyco.dat or .par file.
Note that if you give it a .par file it will STILL use the most recent polyco_new.dat UNLESS you specify -r|--recalc-polycos!
This means it is convenient to ALWAYS give it the .par file name and just use -r when you want to recalc."""
desc=desc+"\n\n"+"\t\t\t ** Emilie Parent -- Modifications, 2022-02: GEOCENTERED or BARYCENTERED events from X-ray satellites.**"

parser=OptionParser(usage=" %prog [options] [FT1_FILENAME] [POLYCO_FILENAME/PAR_FILENAME]",
                                       description=desc)
parser.add_option("-n","--ntoa",type="int",default=10,help="Number of TOAs to produce between TSTART and TSTOP.")
parser.add_option("-b","--nbins",type="int",default=32,help="Number of bins in each profile.")
parser.add_option("-e","--emin",type="float",default=0.0,help="Minimum energy to include.")
parser.add_option("-x","--emax",type="float",default=300000.0,help="Maximum energy to include.")
parser.add_option("-p","--plot",type="string",default=None,help="Set to base name of files to generate likelihood surface plots for each TOA.")
parser.add_option("","--wmin",type="float",default=0.0,help="Minimum weight value to include.")
parser.add_option("","--wmax",type="float",default=1.0,help="Maximum weight value to include.")
parser.add_option("-t","--template",type="string",default=None,help="File name containing template (LCTemplate compatible)")
parser.add_option("-w","--weights",type="string",default=None,help="Specify column of FT1 file to use as weights for TOA computation")
parser.add_option("-r","--recalc-polycos",action="store_true",default=False,help="Recompute polycos using Tempo2")
parser.add_option("-u","--uniform_sigma",action="store_true",default=False,help="Instead of fixed time bins, make a TOA each time the H-test gets to a sufficient sigma.")
parser.add_option("-a","--addphase",action="store_true",default=False,help="Add PULSE_PHASE column to FT1 file using poly cos.")
parser.add_option("","--blind",action="store_true",default=False,help="Force blind search for TOAs rather than tracking.")
parser.add_option("-o","--output",type="string",default=None,help="File for output of .tim file.  Otherwise output to STDOUT.")
parser.add_option("","--calc_photonphase",action="store_true",default=False,help="Calculate photon phases with PINT's photonphase program "\
								"and writes PULSE_PHASE header entries to new file (extension _PHASE.fits).")

## Parse arguments
(options,args) = parser.parse_args()
if len(args) != 2:
    parser.print_help()
    parser.error("Both FT1_FILENAME and POLYCO_FILENAME arguments are required.")

#if options.gauss is None and options.unbinned:
#    print >>sys.stderr, "ERROR! Unbinned likelihood requires a gaussian template."

ft1name    = args[0]
parfile    = args[1]

start_message(ft1name)

timing_model = models.get_model(parfile)


# Call PINT's photonphase 
if options.calc_photonphase or not (has_phase_info(ft1name)):
    extension = '.'+ft1name.split('.')[-1]
    outfile = ft1name.split(extension)[0] + '_PHASE.fits'
    outplot = ft1name.split(extension)[0] + '_photonfold.png'
    def_args = ['--addphase','--absphase','--barytime','--plot', '--plotfile',outplot,'--outfile',outfile,ft1name,parfile]
    print("\n ** Calling PINT's photonphase program to calculate photon phases\n")
    photonphasemain(def_args)
    ft1name = outfile  
    
    print("\nTOAs will now be extracted from file: \n %s\n"%outfile)

timesys, site = get_time_sys(ft1name)
poly = polycos.Polycos()

print(" Reading in data ... " )
data = PhaseData(ft1name,poly,use_weights=(options.weights is not None),we_col_name=options.weights,wmin=options.wmin,emin=options.emin,emax=options.emax)

print(" Computing polycos ... ")
segL = int(np.ceil((data.mjd_stop - data.mjd_start)*24*60./options.ntoa))
poly.generate_polycos(model=timing_model, mjdStart=data.mjd_start-0.5, mjdEnd=data.mjd_stop+0.5, obs=site, segLength=segL, ncoeff=4, obsFreq=0)

if options.addphase: data.write_phase()

#template = LCTemplate(template=options.template)
norms = None
prims = lctemplate.prim_io(template=options.template)		# EP: use pint's "prim_io" and LCTemplate instead of GeoTOA's
if len(prims)>1:
	norms = prims[-1]
	prims = prims[0]
template = lctemplate.LCTemplate(primitives=prims,norms=norms)

if options.weights is not None:
    logl = np.log(1+data.weights*(template(data.ph)-1)).sum()
else: logl = np.log(template(data.ph)).sum()
print('\n Total log likelihood:  %.2f\n'%(logl))


print(" Extracting TOAs ... ")
tg = UnbinnedTOAGenerator(data,poly,template,plot_stem=options.plot,good_ephemeris=(not options.blind))


if options.uniform_sigma:
    binner = UniformLogBinner(options.ntoa,data,template)
else:
    binner = TOABinner(options.ntoa,data)
    
toas,err_toas,tim_strings = tg.get_toas(binner)

if options.output is not None:
    f = open(options.output,'w')
    f.write('\n'.join(tim_strings))
    f.write('\n')
    f.close()

    print('Wrote princeton-format TOAs to %s'%(options.output))

else: print('\n'.join(tim_strings))