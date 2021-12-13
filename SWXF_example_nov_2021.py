# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 17:47:49 2021

SWXF_example_nov_2021

An example file as a first pass at simulating x-ray
standing wave fluorescence from labeled proteins
on a solid support.  

Laurence Lurio, Northern Illinois University
December 13, 2021

@author: lluri
"""
# pip install xraydb
# use MKS units
import xraydb as xdb
import numpy as np
import scipy.constants as scc
from matplotlib import pyplot as plt



# Define some constants use MKS units
nm = 1e-9
cm = .01
micro = 1e-6
pi =  scc.pi
degree = pi/180
N_A = scc.physical_constants["Avogadro constant"][0]
r_0 = scc.physical_constants["classical electron radius"][0]
q_e = scc.physical_constants["elementary charge"][0]
h = scc.physical_constants["Planck constant"][0]
c = scc.physical_constants["speed of light in vacuum"][0]

# Define some material properties
D = 5*nm                # Diameter of gold nanoparticles
V = (1/6)*pi*D**3       # Volume of gold nanoparticles
# Calculate number of gold atoms in nanoparticle
N_gold = V*xdb.atomic_density('Au')*N_A/xdb.atomic_mass('Au')/cm**3
# now assume a density of one nanoparticle for every 100 nm^2 and find the
# surface coverage of gold
print('Number of gold atoms per nanoparticle = {0:7.0f}'.format(N_gold))
Coverage = 1/100**2/nm**2       #number of gold nanoparticles per nm^2
N_per_A = N_gold*Coverage
print('Number of gold atoms per m^2 = {0:7.3g}'.format(N_per_A))
# Next find the fluorescence cross section per unit area on the surface
# Pick an energy just above the fluorescence  L edge of gold
E = 12000       # Energies are assumed to be in eV
# Calculate the cross section per atom
sig = xdb.mu_elam('Au', E, kind='photo')*xdb.atomic_mass('Au')*cm**2/N_A
# Correct for the fluorescence yield
fyield = xdb.xray_edges('Au')['L3'][1]
sig_total_per_meter= sig*N_per_A*cm**2*fyield
print('fractional cross section per meter = {0:7.3g}'.format(sig_total_per_meter))
# Assume an incident beam with 10^10 photons/s
I0 = 1e10
# Assume an incident angle of 0.1 degree
# First calculation, ignore reflection, just calculate fluorescence yield
# assuming beam comes in at a small grazing angle and doesn't reflect
alpha = 0.1*degree
A_beam = (20*micro)**2
phi0 = I0/A_beam
sig_total = sig_total_per_meter*A_beam/alpha
yeild = sig_total*phi0
print('expect total yeild of {0:7.3f} photons/s above critical angle'.format(yeild))
# Now define a function to represent the standing wave assuming specular reflectivity 
# from the surface using Fresnels equations for x-rays
#
# Take the mirror surface to be Silicon (element #14)
#
# First define a function which calculates the x-ray index of refraction for
# the mirror
def n_elem(elem,E): 
    # subroutine to calculate the x-ray index of refraction of 
    # an elemental material
    rho = xdb.atomic_density(elem)
    f = xdb.atomic_number(elem) + xdb.f1_chantler(elem,E)
    Ne = rho*f*N_A*1e6/xdb.atomic_mass(elem)    
    lam = h*c/q_e/E
    n = 1.0 - r_0*Ne*lam**2/2.0/scc.pi
    return n
def n_water(E):
    # subroutine to calculate the x-ray index of refraction of water  
    f1_H = xdb.atomic_number('H') + xdb.f1_chantler('H',E)
    f2_H = xdb.f2_chantler('H',E)
    f1_O = xdb.atomic_number('O') + xdb.f1_chantler('O',E)
    f2_O = xdb.f2_chantler('O',E)
    A_H = xdb.atomic_mass('H') 
    A_O = xdb.atomic_mass('O') 
    f = 2*(f1_H+1j*f2_H)+(f1_O+1j*f2_O)
    A = A_O+2*A_H
    rho = 1.0 # density of water 1 g/cc
    Ne = rho*f*N_A*1e6/A    
    lam = h*c/q_e/E
    n = 1.0 - r_0*Ne*lam**2/2.0/scc.pi
    return n

def swave(alpha,z,E):
    k0 = 2*pi*E*q_e/h/c
    thc_si = np.sqrt(2*(1-n_elem(14,E)))
    thc_water = np.sqrt(2*(1-n_water(E)))
    thc=np.sqrt(thc_si**2-thc_water**2)
    alphap = np.sqrt(alpha**2-thc**2)
    R=(alpha-alphap)/(alpha+alphap)
    Er=R*np.exp(-1j*k0*z*np.sin(alpha))
    E0=np.exp(1j*k0*z*np.sin(alpha))
    return abs(E0+Er)**2

#%% Now make some plots of the standing wave
# Now calculate standing wave
plt.figure(1)
plt.clf()
z=np. linspace(0,100)*nm # setup an array of height positions
alpha = 0.08*degree # incident angle
Etot = swave(alpha,z,E)
plt.plot(z/nm,Etot)
# indicate position of Au nanoparticles 5 nm above the surface
z_au = 5*nm
plt.figure(1)
plt.clf()
plt.plot(z/nm,Etot)
plt.plot(np.array([z_au,z_au])/nm,np.array([0,4]),'--r')
plt.xlabel('height above surface (nm)')
plt.ylabel('$I/I_0$')
plt.title('Standing wave amplitude for {0:2.0f} KeV wave above Si surface at {1:3.2f} deg'.format(
    E/1000,alpha/degree))
# second plot calculate fluorescence yeild vs. angle
alpha = np.linspace(0.01,2,1000)*degree
Iz = swave(alpha,z_au,E)
fyeild = Iz*sig_total_per_meter*A_beam*phi0/(alpha*0+.1*degree)
plt.figure(2)
plt.clf()
plt.plot(alpha/degree,fyeild)