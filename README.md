# High-Energy TOA extraction tool for pulsar timing purposes

Package: HighE-TOA (adapted from GeoTOA by Matthew Kerr & Paul Ray, see below) \
Authors: Emilie Parent <parent@ice.csic.es> \
Version: 1.0 \
Date:	   2022 February 28 

This package is an adaptation of the Geocentric TOA Tools (GeoTOA): <https://fermi.gsfc.nasa.gov/ssc/data/analysis/user/>
developed by Matthew Kerr <matthew.kerr@nrl.navy.mil> and Paul Ray <paul.ray@nrl.navy.mil>. \
  *Package: Geocentric TOA Tools (GeoTOA)* \
  *Authors: Matthew Kerr <matthew.kerr@nrl.navy.mil> and Paul Ray <paul.ray@nrl.navy.mil>* \
  *Version: 2.0* \
  *Date: 2022 January 11* 

***NOTICE***\
This code is provided as contributed code in hopes that
it will be helpful. The authors do not warrant that the code is bug free,
or even correct. E-mail questions and requests are welcomed, but support
can only be provided on a time-available basis.

The main differences between the HighE-TOA and GeoTOA packages are:

* Can read-in barycentered TOAs from other high-energy instruments (whereas GeoTOA
  only takes-in geocentered LAT data)

* Some python modules and functions in GeoTOA have corresponding implementation in 
  PINT: HighE-TOA imports those modules and functions directly from PINT 

* Option to calculate photon phase with PINT's photonphase tool (implemented in both
  itemplate.py and upolyfold.py) and add a PULSE_PHASE column to input data. Using tempo2 on
  raw data (required with GeoTOA) is no longer required with HighE-TOA. 

* Option to specify the number of profile bins for folding (implemented in both
  itemplate.py and upolyfold.py).

* Template profile can be fitted with Gaussians using a modified version of 
  presto's pygaussfit (presto_pygaussfit.py).

* By default, TOAs are written in princeton format, which is readible by both TEMPO and TEMPO2 

----------------------------------------------------------
The following instructions have been adapted from [GeoTOA's](https://fermi.gsfc.nasa.gov/ssc/data/analysis/user/GeoTOA_README.txt): 

**PREREQUISITES:**
* [Fermitools](https://github.com/fermi-lat/Fermitools-conda/wiki) software, version 2.0.8 or later (Python 3 version)
* [Tempo](http://tempo.sourceforge.net/) or [Tempo2](http://www.atnf.csiro.au/research/pulsar/tempo2/) software installed 
* A version of python with numpy, matplotlib and scipy installed.  
* [PINT](https://github.com/nanograv/PINT)  pulsar timing software
* [PRESTO](https://github.com/scottransom/presto) suite of pulsar software ** Only required if fitting template using 'presto' method.

**INSTALLATION:**
* Clone HighE-TOA (or download source) from github. This will make a directory called HighE-TOA with two subdirectories (bin and python)
* Add the bin directory to your PATH
* Add the python directory to your PYTHONPATH, typically like this:\
 ` setenv PYTHONPATH <installdir>/HighE-TOA/python:$PYTHONPATH `

**DESCRIPTION OF CODES:**
* upolyfold.py computes TOAs from unbinned data (a photon, or 'FT1' file) using the maximum likelihood methods described in [Ray et al. 2011, ApJS, 194, 17](https://ui.adsabs.harvard.edu/abs/2011ApJS..194...17R/abstract).
  Will also compute photon phases if input file does not contain 'PULSE_PHASE' column.
* itemplate.py generates a template for use with upolyfold.py. Will also compute photon phases if input file does not contain 'PULSE_PHASE' column.
* presto_pygaussfit.py: interactive profile fitter (multi-component Gaussians). Will write the results to file <filename>prof_16bins_template.gaussian,
  which can then be used as template when calling upolyfold.py 

**REFERENCES AND CITATIONS:** \
 This package use routines from the Geocentric TOA Tools ([GeoTOA](https://fermi.gsfc.nasa.gov/ssc/data/analysis/user)) 
 by Matthew Kerr and Paul Ray. The unbinned maximum likelihood methods are described in [Ray et al. 2011, ApJS, 194, 17](https://ui.adsabs.harvard.edu/abs/2011ApJS..194...17R/abstract).\


**USAGE:**

* Make sure you have a photon event file corresponding to the energy
  and angle cuts that maximize the S/N for your pulsar (GeoTOA: event file 
  extracted from the FSSC archive, and the spacecraft file that covers the 
  same time range -- no longer necessary with HighE-TOA.)

* You also need a .par file with a timing model for your pulsar that
  works at least well enough that you can fold the pulsar without
  losing phase count.

* (GeoTOA: Use tempo2 -gr fermi to add a PULSE_PHASE column to your photon
  file. To use Tempo2 to do this you MUST use the raw (NOT barycentered
  or geocentered) photon event file.)

* Use itemplate.py to generate a template profile for your
  timing.  Just run: \
  ` itemplate.py --nbins <N> my_photon_file.fits> ` \
  and tell it you want a gaussian template (with either 'gauss' or 'presto' method). 
  Then draw in some gaussian components by clicking and dragging your 
  mouse across the FWHM of each visible component.  Close the window and 
  itemplate will do an unbinned fit to the model. If using 'gauss' method: give 
  it a file name to write out your template (for example "template.gauss")
  (automatic if using 'presto' method). \
  If input file does not have a 'PULSE_PHASE' column, itemplate.py will first offer
  to compute photon phases (will ask user for .par file for phase prediction) 
  and will write into new file <my_photon_file>_PHASE.fits. Will then proceed with 
  template generation. 
  
* (GeoTOA: Now, to generate TOAs, you need a photon file with geocentered times.
  Generate this using "gtbary -tcorrect=GEO". Give the output file a
  name that indicated it contains geocentered data, so you don't get
  confused!) -- No longer necesary with HighE-TOA. Go straight to next step. 

* Type upolyfold.py -h and read about the options!

* Then, to create TOAs using upolyfold.py, use a command like this:\
` upolyfold.py -r -b <N> -n 32 -t template.3gauss J1231-1411_1.0Deg_239557517_361106549MET_350_300000MeV_z100_GEO.fits J1231-1411.par `

This will create 32 TOAs from your photon file using the template called
"template.3gauss".  The folding will be done with polycos generated
from "J1231-1411.par".  The polycos get stored in polyco_new.dat. If
you want to reuse that file without recomputing it, just leave the
"-r" off the command line.  This CAN BE UNSAFE because if you have
changed the par file in the meantime, the polycos will be old.  Be
careful!

You can also add "-o outfile.tim" to write the TOAs to an output file.
Each TOA is comments with the number of photons that went into that
TOA ("-np #") and also the probability that the observed photon
distribution in that TOA's data came from a uniform (unpulsed)
sources.  These can be helpful diagnostics of whether a TOA is good.
