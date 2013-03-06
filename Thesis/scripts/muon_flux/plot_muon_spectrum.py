#!/bin/env python

import ROOT
import numpy
import array
import scipy.integrate
from math import pi, sin, cos, acos, pow, exp, log

def surface_flux(theta, E):
    if E < 0:
        return 0
    else:
        return 0.14*pow(E, -2.7)*(1/(1+1.1*E*cos(theta)/115)
                                  +0.054/(1+1.1*E*cos(theta)/850))

def get_a(E):
    # GeV cm^2/g
    if E <= 10:
        return 1e-3*2.17
    elif E > 10 and E <= 100:
        return 1e-3*(2.17 + (2.44-2.17)*(E-10)/(100-10))
    elif E > 100 and E <= 1000:
        return 1e-3*(2.44 + (2.68-2.44)*(E-10)/(1000-100))
    elif E > 1000 and E <= 10000:
        return 1e-3*(2.68 + (2.93-2.68)*(E-10)/(10000-1000))
    else:
        return 1e-3*2.93

def get_b(E):
    # cm^2/g
    if E <= 10:
        return 1e-6*1.90
    elif E > 10 and E <= 100:
        return 1e-6*(1.90 + (3.04-1.90)*(E-10)/(100-10))
    elif E > 100 and E <= 1000:
         return 1e-6*(3.04 + (3.92-3.04)*(E-100)/(1000-100))
    elif E > 1000 and E <= 10000:
        return 1e-6*(3.92 + (4.35-3.92)*(E-1000)/(10000-1000))
    else:
        return 1e-6*(4.35)

def flux(theta, h, E):
    # Slant Depth
    X = h*100/cos(theta)  # g/cm^2 (1 m.w.e. = 1e2 g/cm^2)

    a = get_a(E)
    b = get_b(E)

    eps = a/b
    E_0 = E*exp(b*X) + eps*(exp(b*X)-1)

    return surface_flux(theta, E_0)*exp(b*X)

def flux_sin_theta(theta, h, E):
    return sin(theta)*flux(theta, h, E)

def integrated_flux(E):
    return 2*pi*scipy.integrate.quad(flux_sin_theta, 0, pi/2, (1526, E))[0]

def meihime_flux(E):
    h = 1526*100 # g/cm^2 (1 m.w.e. = 1e2 g/cm^2)

    a = get_a(E)
    b = get_b(E)
    eps = a/b

    g = 3.77

    return exp(-b*h*(g-1))*pow(E+eps*(1-exp(-b*h)), -g)

def miyake_ang(theta):
    h = 1585*100 # g/cm^2 (1 m.w.e. = 1e2 g/cm^2)

    return sin(theta)*pow(cos(theta), 1.53)*exp(-8e-6*h*(1/cos(theta)-1))

def meihime_ang(theta):
    h = 1585*1e-3 # 1e5g/cm^2

    X = h

    l1 = 0.45
    l2 = 0.87
    I1 = 8.6e-6
    I2 = 0.44e-6

    return sin(theta)*((I1*exp(-X/(l1*cos(theta)))+I2*exp(-X/(l2*cos(theta))))/cos(theta))
    

if __name__ == "__main__":

    Es = numpy.logspace(0, 4, 1000, False)
    
    E = array.array('d', Es)
    norm = integrated_flux(Es[0])
    N_g = array.array('d', [integrated_flux(x)/norm for x in Es])
    norm = meihime_flux(Es[0])
    N_mh = array.array('d', [meihime_flux(x)/norm for x in Es])
    norm = flux(0, 1500, Es[0])
    N_00 = array.array('d', [flux(0, 1500, x)/norm for x in Es])
    N_30 = array.array('d', [flux(pi/6, 1500, x)/norm for x in Es])
    N_60 = array.array('d', [flux(pi/3, 1500, x)/norm for x in Es])

    #print("Maximum: %e."%(max(N)))
    #print("Maximum at E=%e."%(Es[N.index(max(N))]))

    c_E = ROOT.TCanvas("c_E")
    g_E_g = ROOT.TGraph(len(E), E, N_g)
    g_E_g.SetTitle("Integrated Flux;E (GeV);normalized flux (GeV^{-1}")
    g_E_mh = ROOT.TGraph(len(E), E, N_mh)
    g_E_mh.SetLineColor(ROOT.kRed)
    g_E_g.Draw("AL")
    g_E_mh.Draw("L")
    c_E.SetLogy()
    c_E.SetLogx()
    c_E.Update()

    c_ang = ROOT.TCanvas("c_ang")
    g_ang = ROOT.TMultiGraph("g_ang",
                             "Angular Flux;E (GeV); normalized flux (GeV^{-1})")

    g_00 = ROOT.TGraph(len(E), E, N_00)
    g_30 = ROOT.TGraph(len(E), E, N_30)
    g_60 = ROOT.TGraph(len(E), E, N_60)
    g_30.SetLineColor(ROOT.kBlue)
    g_60.SetLineColor(ROOT.kRed)
    g_ang.Add(g_00, "l")
    g_ang.Add(g_30, "l")
    g_ang.Add(g_60, "l")
    l_ang = ROOT.TLegend(0.11, 0.11, 0.26, 0.26)
    l_ang.SetFillColor(ROOT.kWhite)
    l_ang.AddEntry(g_00, "#theta = 0", "l")
    l_ang.AddEntry(g_00, "#theta = #pi/6", "l")
    l_ang.AddEntry(g_60, "#theta = #pi/3", "l")
    g_ang.Draw("AL")
    l_ang.Draw()
    c_ang.SetLogy()
    c_ang.SetLogx()
    c_ang.Update()

    theta = array.array('d', numpy.linspace(0, pi/2, 1000, False))
    norm = scipy.integrate.quad(miyake_ang, 0, pi/2)[0]
    norm_mi = norm
    #norm = miyake_ang(acos(1))
    N_ang_mi = array.array('d', [miyake_ang(x)/norm for x in theta])
    norm = scipy.integrate.quad(meihime_ang, 0, pi/2)[0]
    #norm = meihime_ang(acos(1))
    N_ang_mh = array.array('d', [meihime_ang(x)/norm for x in theta])

    print("Maximum: %e."%(max([miyake_ang(x) for x in theta])))
    #print("Maximum at E=%e."%(Es[N.index(max(N))]))

    c_theta = ROOT.TCanvas("c_theta")
    g_mh =  ROOT.TGraph(len(theta), theta, N_ang_mh)
    g_mh.SetTitle("Angular Distribution;cos(#theta);normalized flux")
    g_mi =  ROOT.TGraph(len(theta), theta, N_ang_mi)
    g_mh.SetLineColor(ROOT.kRed)
    g_mh.Draw("AL")
    g_mi.Draw("L")
    c_theta.Update()
    
    dummy = raw_input("Press Enter...")
