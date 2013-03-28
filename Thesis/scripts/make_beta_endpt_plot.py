import ROOT
import numpy
from math import sqrt, pi, exp

alpha = 1.0/137
m_e = 511e3 # eV
Q = 18574 # eV

def p(E):
    # relativistic momentum for KE E
    return sqrt((E + m_e)**2 - m_e**2)

def beta(E):
    # relativistic beta for KE E
    return p(E)/(E + m_e)

def A(E, Z):
    # for use in Fermi function
    return 2*pi*alpha*(Z+1)/beta(E)

def F(E, Z):
    # Fermi function
    x = A(E, Z)
    return x / (1 - exp(-x))

def R(E):
     # spectrum for tritium dR/dE
     Z = 1
     return F(E, Z) * p(E) * (E + m_e) * (Q - E)

def R_mass(E, m_nu):
    # spectrum considering neutrino mass
    try:
        return R(E) * sqrt((Q-E)**2 - m_nu**2)
    except ValueError:
        return 0

scale = 1.0/R_mass(0.9999*Q, 0)

mg = ROOT.TMultiGraph("mg_endpt", "")
l = ROOT.TLegend(0.78, 0.78, 0.98, 0.98)
l.SetFillColor(ROOT.kWhite)
graphs = []
colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue]

for (m_nu, color) in zip([0, 1, 0.25], colors):
    g = ROOT.TGraph(0)

    for E in numpy.linspace(0.9999, 1, 1000, False):
        n = g.GetN()
        g.Set(n+1)
        r = scale*R_mass(E*Q, m_nu)
        #if r == 0:
        #    r = 1e-6
        g.SetPoint(n, E, r)

    g.SetLineColor(color)
    g.SetTitle("m_{#nu} = %0.2f eV"%(m_nu))
    graphs.append(g)
    l.AddEntry(g, "m_{#nu} = %0.2f eV"%(m_nu), "l")
    mg.Add(g, "l")

c = ROOT.TCanvas("c")
c.SetRightMargin(0.01)
c.SetTopMargin(0.01)

mg.Draw("AL")
mg.GetXaxis().SetTitle("E/Q")
mg.GetYaxis().SetTitle("dN/dE (relative)")
l.Draw()

#c.SetLogy()

dummy = raw_input("Press Enter...")
