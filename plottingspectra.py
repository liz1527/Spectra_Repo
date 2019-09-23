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
        'size'   : 14}

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

DR1_IDs = ['77042', '66065', '43445','44239','48174','51911',]
DR11_IDs = ['252446', '219656', '148305', '148895','161351','173520',]
zs = [3.2176, 0.4898, 0.4815,0.916,1.4102,2.9982]


def drawspectrumUDSz(DR1_ID, DR11_ID, z):
    
    spec = read_fits.read_fits_spectrum1d('variable-spectra/UDSz_DR1_'+DR1_ID+'_DR11_'+DR11_ID+'.fits', dispersion_unit='angstrom')
    plt.figure(figsize=[17,4])
    plt.plot(spec.wavelength, spec.flux)
    axes = plt.gca()
    ylims = axes.get_ylim()
    
    for key in linedict:
        val = linedict[key]
        lineval = val*(z+1)*u.AA
        if lineval < max(spec.wavelength) and lineval > min(spec.wavelength):
            plt.vlines(lineval.value, ymin=ylims[0], ymax=ylims[1], linestyles='dashed')
            if key == 'SiII' or key =='LyA' or key == 'LyB':
                plt.text(lineval.value-140, ylims[0]+0.2e-18, key)
            elif key == 'SiIV_1393+OIV1402':
                plt.text(lineval.value-400, ylims[0]+0.2e-18, key)
            elif key == 'Balmer_break':
                plt.text(lineval.value-300, ylims[1]+0.2e-18, key)
            elif key == 'OIII_doublet-1/3':
                plt.text(lineval.value-200, ylims[1]+0.2e-18, 'OIII_doublet')
            elif key == 'OIII_doublet-1':
                continue
            else:
                plt.text(lineval.value+10, ylims[0]+0.2e-18, key)
        
    plt.ylim(ylims)
    plt.ylabel('Flux')
    unit = u.AA
    plt.xlabel('Wavelength, '+unit.to_string('latex'))
    plt.tight_layout()
    plt.savefig(DR11_ID+'spectra.pdf')
    
def drawspectrum3DHST(DR8_ID, DR11_ID, z):
    
    spec = fits.getdata('variable-spectra/3DHST_DR8_'+DR8_ID+'_DR11_'+DR11_ID+'.fits')
    plt.figure(figsize=[17,4])
    plt.plot(spec.wave, spec.flux)
    axes = plt.gca()
    ylims = axes.get_ylim()
    
    for key in linedict:
        val = linedict[key]
        lineval = val*(z+1)
        if lineval < max(spec.wave) and lineval > min(spec.wave):
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
    plt.tight_layout()
#    plt.savefig(DR11_ID+'spectra.png')
    
#drawspectrumUDSz(DR1_IDs[0], DR11_IDs[0], zs[0])
#drawspectrumUDSz(DR1_IDs[1], DR11_IDs[1], zs[1])
#drawspectrumUDSz(DR1_IDs[2], DR11_IDs[2], zs[2])
#drawspectrumUDSz(DR1_IDs[3], DR11_IDs[3], zs[3])
#drawspectrumUDSz(DR1_IDs[4], DR11_IDs[4], zs[4])
#drawspectrumUDSz(DR1_IDs[5], DR11_IDs[5], zs[5])
#drawspectrumUDSz('47304', '159009', 1.7591)
#drawspectrum3DHST('74491', '118344', 0.84688)
#drawspectrum3DHST('59946', '94947', 1.5434)