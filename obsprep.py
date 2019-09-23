#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 10:14:00 2018

@author: ppxee
"""

import matplotlib.pyplot as plt #for plotting
import matplotlib.patches as patches #for square
import numpy as np
from astropy.table import Table
plt.close('all')

plt.figure()
data = Table.read('variablesforfors2maglim23.fits')
plt.scatter(data['RA'], data['DEC'], c=data['RMAG_20'])
plt.gca().invert_xaxis()
plt.gca().add_patch(patches.Rectangle((34.4, -5.1), 6.83/60, 6.83/60, fill=False))      # (x,y)# width # height
plt.gca().axis('equal')
plt.xlabel('RA')
plt.ylabel('Dec')
plt.colorbar()

#plt.figure()
#plt.scatter(data['RA'], data['DEC'], c=data['RMAG_20'])
#plt.gca().invert_xaxis()
#plt.gca().add_patch(patches.Rectangle((34.4, -5.35), 27.4/60, 27.4/60, fill=False))      # (x,y)# width # height
#plt.gca().axis('equal')
#plt.colorbar()
#
#
#plt.figure()
data2 = Table.read('variablesforfors2thresh3.5noxraymaglim23.fits')
#plt.scatter(data2['RA'], data2['DEC'], c=data2['RMAG_20'])
#plt.gca().invert_xaxis()
#plt.gca().add_patch(patches.Rectangle((34.4, -5.1), 6.83/60, 6.83/60, fill=False))      # (x,y)# width # height
#plt.gca().axis('equal')
#plt.colorbar()
#
plt.figure()
plt.scatter(data2['RA'], data2['DEC'], c=data2['RMAG_20'])
plt.gca().invert_xaxis()
plt.gca().add_patch(patches.Rectangle((34.4, -5.35), 27.4/60, 27.4/60, fill=False))      # (x,y)# width # height
plt.gca().axis('equal')
plt.xlabel('RA')
plt.ylabel('Dec')
plt.colorbar()

plt.figure()
data2 = Table.read('variablesforfors2thresh3.5noxraymaglim23.fits')
data3 = Table.read('variablesforfors2thresh3.5maglim23.fits')
plt.scatter(data3['RA'], data3['DEC'], c='r')
plt.scatter(data2['RA'], data2['DEC'], c='k')
plt.scatter(data['RA'], data['DEC'], c='b')
plt.gca().axis('equal')
plt.gca().invert_xaxis()
plt.xlabel('RA')
plt.ylabel('Dec')
plt.gca().add_patch(patches.Rectangle((34.4, -5.1), 6.83/60, 6.83/60, fill=False))      # (x,y)# width # height
#plt.colorbar()