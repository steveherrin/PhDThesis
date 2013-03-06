#!/bin/env python

import scipy.optimize
from math import exp

# using 1526 MC
vertical_flux = 4.24e-07 # Hz/(cm^2 sr)
vertical_flux_err = 5e-9

# Esch paper
#vertical_flux = 3.1e-07 # Hz/(cm^2 sr)
#vertical_flux_err = 0.06e-7

def flux(h):

    A1 = -11.22
    A2 = -0.00262
    A3 = -14.1
    A4 = -0.001213
    A5 = 2.17e-13

    return exp(A1 + A2*h) + exp(A3 + A4*h) + A5

def flux_meihime(h):

    A1 = 8.6e-6
    l1 = 450
    A2 = 0.44e-6
    l2 = 870

    return A1*exp(-h/l1) + A2*exp(-h/l2)


def f(h, vflux):
    
    return flux(h) - vflux

h = scipy.optimize.brentq(f, 0, 1e6, (vertical_flux))

h_pos = scipy.optimize.brentq(f, 0, 1e6,
                              (vertical_flux - vertical_flux_err))

h_neg = scipy.optimize.brentq(f, 0, 1e6,
                              (vertical_flux + vertical_flux_err))

print("A flux of (%.2e +/- %.2e) Hz/(cm^2 sr)"
      %(vertical_flux, vertical_flux_err)
      + " corresponds to a depth of (%0.1f +%0.1f -%0.1f) m.w.e."
      %(h, h_pos-h, h-h_neg))
