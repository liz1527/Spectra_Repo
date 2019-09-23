#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 10:25:50 2018

@author: ppxee
"""

from specutils.io import read_fits
from astropy import units as u
import matplotlib.pyplot as plt #for plotting
import numpy as np
from astropy.io import fits
plt.close('all')

font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 18}

plt.rc('font', **font)
### Set up a dictionary of important lines ###

linedict = {'CIV': 1549.00,
            'SiII': 1526.70,
            'FeII': 1608.50,
            'CIII': 1909.00,
            'MgII': 2799.00,
            'OII': 3727.50,
            'OIII_doublet-1': 5006.80,
            'OIII_doublet-1/3': 4958.90,
            'MgI': 5175.40,
            'LyA': 1215.70,
            'Hbeta': 4861.30,
            'Hgamma': 4340.40,
            'Hdelta': 4101.70,
            'Balmer_break': 4000,
            'NV': 1240.00,
            'SiIV_1393+OIV1402': 1397.00,
            'LyB': 1025.60,
            'OVI': 1033.01}

spec_varys = fits.open('no06_variables_chi30_2arcsec_spec_DR11.fits')[1].data
DR11_ID = spec_varys['DR11-ID']
spec_varys = spec_varys[DR11_ID==93021]#93021, 206729, 217401
cat_ID = spec_varys['cat-ID']
zspec = spec_varys['z-spec']
DR11_ID = spec_varys['DR11-ID']
xray = spec_varys['X-ray']

def drawspectrum(cat_ID, DR11_ID, z, Xray):
    spec = fits.open('UDS-spectra-Feb-19/UDS-'+str(cat_ID)+'-1Dspec.fits')[1].data
    plt.figure(figsize=[20,4])
    if Xray == True:    
        plt.plot(spec['wavelength'], spec['flux'],'r')
    else:
        plt.plot(spec['wavelength'], spec['flux'],'b')
#    
#    plt.plot(spec['wavelength'], spec['flux'],'r')
#    plt.plot(spec['wavelength'], spec['flux'],'b')

    axes = plt.gca()
    plt.ylim(ymax=1.7e-18, ymin=-0.6e-18) #93021
    plt.xlim(xmax=9500, xmin=4800) #93021
#    plt.ylim(ymax=2.6e-18, ymin=-1.2e-18) #206729
#    plt.ylim(ymax=3e-18, ymin=-0.6e-18) #217401
    ylims = axes.get_ylim()
#    xlims = axes.get_xlim()
    
    for key in linedict:
        val = linedict[key]
        lineval = val*(z+1)
        if lineval < max(spec['wavelength']) and lineval > min(spec['wavelength']):
            plt.vlines(lineval, ymin=ylims[0], ymax=ylims[1], linestyles='dashed')
            if key == 'SiII' or key =='LyA' or key == 'LyB':
                plt.text(lineval-140, ylims[0]+0.2e-18, key)
            elif key == 'SiIV_1393+OIV1402':
                plt.text(lineval-400, ylims[0]+0.2e-18, key)
            elif key == 'Balmer_break':
                plt.text(lineval-300, ylims[1]+0.2e-18, key)
            elif key == 'OIII_doublet-1/3':
                plt.text(lineval-200, ylims[1]+0.2e-18, 'OIII_doublet')
            elif key == 'OIII_doublet-1':
                continue
            else:
                plt.text(lineval+10, ylims[0]+0.2e-18, key)
        
    plt.ylim(ylims)
    plt.ylabel('Flux')
    unit = u.AA
    plt.xlabel('Wavelength, '+unit.to_string('latex'))
    plt.title('z = '+str(z))
    plt.tight_layout()
#    plt.savefig('paper_spectra/'+str(DR11_ID)+'spectra.png')
#    plt.savefig('paper_spectra/'+str(DR11_ID)+'spectra.pdf')
#    plt.close('all')
#    
for ID, DR11ID, z, Xray in zip(cat_ID, DR11_ID, zspec,xray):
    drawspectrum(ID, DR11ID, z, Xray)
    
#drawspectrum(cat_ID[-1], DR11_ID[-1], zspec[-1],xray[-1])






