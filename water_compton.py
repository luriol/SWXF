# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 12:18:45 2021
Last updated May 13, 2021
Code to calculate the diffuse scattering from water
and use this to find normalization constant for scattering

@author: Laurence Lurio
"""

import numpy as np
import xraydb
from scipy import constants as scc

# use xraydb to find f0 (Warren calls this fe)
# use O2- as the ion, (2- means two extra electrons)
# thus we ignore the H+ ion
# since it has no electrons
# define scattering vector
E = 13000 # 13000 eV energy 
lam = scc.h*scc.c/scc.e/E/scc.angstrom
print('X-ray wavelength at E = {0:5.0f} eV = {1:5.2f} Angstrom'.format(
    E,lam))
theta = 30*np.pi/180. # scattering angle in radians
# xraydb uses the parameter s to define the 
# scattering wavevector exchange
s = np.sin(theta)/lam
Ne = xraydb.atomic_number('O')+2 # doubly charged O
iM = Ne - xraydb.f0('O2-',s)[0]
print('iM = {0:5.2g}'.format(iM))

# compare with cross section for gold atom
#f'' for gold atom
fpp = xraydb.f2_chantler('Au',E)
print('fpp Au = {0:4.1f}'.format(fpp))
r0 = scc.physical_constants['classical electron radius'][0]
r0 /= scc.angstrom
rat =2*lam*fpp/r0/iM/4/np.pi
print('fluorescence to compton ratio {0:7.3g}'.format(rat))
