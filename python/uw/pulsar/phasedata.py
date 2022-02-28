"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pulsar/phasedata.py,v 1.2 2011/07/21 13:48:49 paulr Exp $

Handle loading of FT1 file and phase folding with polycos.

Authors: Paul S. Ray <paul.ray@nrl.navy.mil>
         Matthew Kerr <matthew.kerr@gmail.com>
"""

import numpy as np
#import pyfits
from astropy.io import fits as pyfits
import sys
from uw.utilities import keyword_options
#from .timeman import METConverter
from astropy import units as u
from astropy import constants as cst
from pint import fits_utils
import psr_utils 

class PhaseData(object):
    """ Object to manage loading data. """

    defaults = (
        ('emin',0,'(MeV) minimum energy to accept'),
        ('emax',99999999,'(MeV) maximum energy to accept'),
        ('we_col_name','WEIGHT','name for column with weight data in FT1 file'),
        ('use_weights',False,'load weights if available'),
        ('wmin',0,'minimum weight to use'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self,ft1file,polyco,**kwargs):
        keyword_options.process(self,kwargs)
        self.ft1file = ft1file
        self.polyco = polyco
        self.process_ft1()

    #def get_freq_info(self):
    #    colnms = self.ft1evt.columns.names
 
    def process_ft1(self):
        #mc = METConverter(self.ft1file)
        
        #self.mjd_start = mc.MJDSTART; self.mjd_stop = mc.MJDSTOP
        f    = pyfits.open(self.ft1file)
        ft1hdr  = f['EVENTS'].header
        ft1evt  = f['EVENTS']
        self.ft1evt  = ft1evt
        self.tel = ft1hdr['TELESCOP']
        srcname = ft1hdr['OBJECT'].replace(' ','').replace('.0-','-').replace('_','')
        self.srcname = srcname[:13]
        #print(self.srcname)
        if 'MJDREF' not in ft1hdr.keys():
            self.mjdref = ft1hdr['MJDREFI']+ft1hdr['MJDREFF']
        else:
            self.mjdref = ft1hdr['MJDREF']
        self.mjd_start = self.mjdref + ft1hdr['TSTART']/86400. + ft1hdr['TIMEZERO']
        self.mjd_stop  = self.mjdref + ft1hdr['TSTOP']/86400. + ft1hdr['TIMEZERO']
        ens  = np.asarray(f['EVENTS'].data.field('PI'))
        self.centre_ens = np.median(ens)
        mask = (ens >= self.emin) & (ens < self.emax)
        if self.use_weights:
            weights = np.asarray(f['EVENTS'].data.field(self.we_col_name))
            mask = mask & (weights > self.wmin)
            self.weights = weights[mask]
        else: self.weights = None
        #mets = np.asarray(f['EVENTS'].data.field('TIME'))[mask]
        #self.mjds = mc(mets)
        if ft1hdr['TIMESYS'] == 'TDB':
            self.BARY = True
            self.obs = '@'
        else:
            self.BARY = False
            self.obs = 'GEO'
        mjds = fits_utils.read_fits_event_mjds(ft1evt)
        self.mjds = np.asarray(mjds)[mask]
        # self.ph = self.polyco.vec_evalphase(self.mjds)
        self.ph = np.asarray(f['EVENTS'].data.field('PULSE_PHASE'))
        #if self.tel != 'NuSTAR':
        try:
            self.get_freq_info()
        #else:
        except:
            self.photon_freq = 0
        f.close()

        #print("Cuts left %d out of %d events." % (mask.sum(), len(mask)), file=sys.stderr)

    def get_freq_info(self):
        colnms = np.asarray(self.ft1evt.columns.names)
        idx = np.where(colnms=='PI')[0][0]
        e_unit =  self.ft1evt.columns[idx].unit
        if e_unit == 'eV':
            photon_e = self.centre_ens*u.eV
        elif e_unit == 'MeV':
            photon_e = self.centre_ens*u.MeV	
        # now convert photon energy to frequency, in MHz
        photon_freq = photon_e/cst.h
        photon_freq = photon_freq.to('MHz')
        self.photon_freq = photon_freq.value

    #def write_phase(self,col_name='PULSE_PHASE'):
    #    f = pyfits.open(self.ft1file)
    #    mc = METConverter(self.ft1file)
    #    mjds = mc(np.asarray(f['EVENTS'].data.field('TIME')))
    #    ph = self.polyco.vec_evalphase(mjds)
        
    #    try:
    #        # clobber old values if there
    #        f['EVENTS'].data.field(col_name)[:] = ph
    #        print('Clobbered old %s column.'%(col_name))
    #    except KeyError:
    #        c    = pyfits.Column(name=col_name,format='D',array=ph)
    #        cols = f['EVENTS'].columns
    #        newcols = cols.add_col(c)
    #        t = pyfits.new_table(newcols,header=f['EVENTS'].header)
    #        t.name = 'EVENTS'
    #        f[1] = t
    #    f.writeto(self.ft1file,clobber=True)

    def toa_data(self,mjd_start,mjd_stop):
        mask = (self.mjds >= mjd_start)&(self.mjds < mjd_stop)
        if self.weights is None: return self.ph[mask],None
        return self.ph[mask],self.weights[mask]
