# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 21:24:16 2022

@author: lluri
"""

from spec_utils import readscan
from matplotlib import pyplot as plt
import numpy as np
specfile = 'align_02222022'
scans = np.array([259, 260, 261, 262])
x = np.array([])
y = np.array([])
for scanno in scans:
    data, scan_info = readscan(specfile,scanno)
    tx = data['Theta'].values
    ty = data['Detector'].values/data['Seconds'].values
    x = np.append(x,tx)
    y = np.append(y,ty)
# plot the combined data
plt.figure('combined plot')
plt.clf()
plt.plot(x,y,'ks')
plt.yscale('log')
plt.xlabel('theta [deg]')
plt.ylabel('counts/second')
plt.title(scan_info['title'])
# simulate reflectivity
def spec_reful(I0,th):
    # Function to approximate reflectivity
    hc = 0.10 # critical angle in degrees
    thf = 3.8  # footprint angle in degrees
    thc = 0.1
    att_water = 0.607
    thp = np.sqrt(th**2-thc**2+0j) 
    I = np.abs((th-thp)/(th+thp))**2
    I *= I0
    I *= th/thf
    I *=att_water
    return I
def rough_layer(th,sig,a,d):
    # function to calcuate effect of an overlayer of thickness d and
    # relative density a
    # both surfaces assumed to have the same roughness sig
    I = th*0+1.0
    lam = 12.398/17.4
    Q = 4*np.pi*np.sin(th/57.3)/lam 
    rough = np.exp(-Q**2*sig**2)
    I *= rough
    lay = (1-a)+a*np.exp(1j*Q*d)
    lay = lay*np.conj(lay)
    I *= np.abs(lay)
    return I
plt.plot(x,spec_reful(1e6,x),'-g',label='Ideal reflectivity')
sig = 4
a = .5
d = 15
plt.plot(x,spec_reful(1e6,x)*rough_layer(x,sig,a,d),'--r',label='with rough hydrocarbon layer')
plt.legend()