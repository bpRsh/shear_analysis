#!python
# -*- coding: utf-8 -*-#
"""
Calculate and plot reduced shear vs other quantities.

@author: Bhishan Poudel

@date:  Feb 2, 2018

"""
# Imports
from __future__ import print_function, division,with_statement,unicode_literals,absolute_import
import sys
import matplotlib.pyplot as plt
import scipy
import numpy as np


# Physics
c         = 3e5    # Speed of light km/h
Hubble    = 67.80  # Hubble constant (km/h)/Mpc
Omega_m   = 0.315  # Current mass density parameter
Omega_d   = 0.685  # Effective density of dark energy

sigma     = 1000   # Lens Velocity dispersion km/h (lens.txt)
z_source  = 0.7    # Source galaxy redshift
z_lens    = 0.3    # Lens redshift
pix_scale = 0.2    # final pixel scale (arcsecond per pixel)


def E(z,Omega_m,Omega_vac): # E(z) as it appears in Hubble parameter evolution
    return np.sqrt(np.reciprocal((Omega_m * (1+z)**3.0)+ Omega_d))


def comoving_dist(E,z1,z2): # Comoving distance
    integral = scipy.integrate.quad(E,z1,z2,args=(Omega_m,Omega_d))
    return c * integral[0] / Hubble

# Angular diameter distance of source from observer
D_s  = comoving_dist(E,0,z_source) / (1+z_source)

# angular diameter distance of souce from lens   
D_ds = comoving_dist(E,z_lens,z_source) / (1+z_source)

# kappa_constant
kappa_constant = 206264.8062471 * (2*np.pi * (sigma*sigma) / (c*c) ) * (D_ds/D_s) / pix_scale



# ellipticity
color_ellip = 'color_galshear_ellip.dat'
mono_ellip =  'mono_galshear_ellip.dat'

# shear
color_shear = 'color_galshear_shear.dat'
mono_shear = 'mono_galshear_shear.dat'

# read values from data files
r    = np.genfromtxt(color_ellip, usecols=(1),unpack=True)  # radius
et_c = np.genfromtxt(color_ellip, usecols=(3),unpack=True)  # et for color
et_m = np.genfromtxt(mono_ellip,  usecols=(3),unpack=True)  # et for mono
gt_c = np.genfromtxt(color_shear, usecols=(3),unpack=True)  # gt for color
gt_m = np.genfromtxt(mono_shear,  usecols=(3),unpack=True)  # gt for mono

# ellipticity and shear ratios
erat = et_c / et_m
grat = gt_c / gt_m

# kappa
kappa = kappa_constant / r

# reduced shear
g = kappa / ( 1- kappa)

# difference and ratio
e_diff = et_c - et_m
e_rat  = et_c / et_m

g_diff = gt_c - gt_m
g_rat  = gt_c / gt_m

# et and em
et_gt = ['et_c', 'et_m', 'e_diff','e_rat', 'gt_c', 'gt_m', 'g_diff', 'g_rat']

# plot
for y in et_gt:
    fig, ax = plt.subplots(figsize=(12,8))
    ax.plot(g,eval(y),ls='-',label=y, lw=1)
    plt.xlabel('g')
    plt.ylabel(y)
    plt.title('Plot of {} vs reduced shear'.format(y))
    plt.savefig('plots/{}_vs_g.pdf'.format(y))

    plt.legend(loc='upper left')
    plt.tight_layout()
    plt.close()
