#!/bin/env python

import ROOT
from math import pi, sqrt, sin, acos, atan, atan2, sin, cos, tan
import array
import numpy

#ROOT.gStyle.SetNumberContours(10)

m_x = ROOT.RooRealVar("#mu_{x}", "#mu_{x}", 1)
m_y = ROOT.RooRealVar("#mu_{y}", "#mu_{y}", 1)

s_x = ROOT.RooRealVar("#sigma_{x}", "#sigma_{x}", 0.07)
s_y = ROOT.RooRealVar("#sigma_{y}", "#sigma_{y}", 0.035)

n_s = 6
E_lo = min(m_x.getVal() - n_s*s_x.getVal(), m_y.getVal() - n_s*s_y.getVal())
E_hi = max(m_x.getVal() + n_s*s_x.getVal(), m_y.getVal() + n_s*s_y.getVal())

E_x = ROOT.RooRealVar("E_x", "Scintillation Energy", E_lo, E_hi, "E/E_{0}")
E_y = ROOT.RooRealVar("E_y", "Ionization Energy", E_lo, E_hi, "E/E_{0}")

rho = ROOT.RooRealVar("#rho", "#rho", -0.80)
theta = ROOT.RooRealVar("#theta", "#theta", acos(s_x.getVal()*(s_x.getVal()-rho.getVal()*s_y.getVal())
                                                 /sqrt(s_x.getVal()**4+s_y.getVal()**4
                                                       -2*rho.getVal()*s_x.getVal()**3*s_y.getVal()
                                                       -2*rho.getVal()*s_x.getVal()*s_y.getVal()**3
                                                       +2*(rho.getVal()*s_x.getVal()*s_y.getVal())**2)))

print("I'm using       theta = %0.2f"%(theta.getVal()*180/pi))
print("Conti suggests  theta = %0.2f"%(0.5*atan(2*rho.getVal()*s_y.getVal()*s_x.getVal()/(s_y.getVal()**2-s_x.getVal()**2))*180/pi))
print("Aprile suggests theta = %0.2f"%(atan2(s_y.getVal(), s_x.getVal())*180/pi))

m_u = ROOT.RooFormulaVar("#mu_{u}", "#mu_{u}", "(@1*sin(@0) + @2*cos(@0))/(sin(@0)+cos(@0))",
                         ROOT.RooArgList(theta, m_x, m_y))
m_v = ROOT.RooFormulaVar("#mu_{v}", "#mu_{v}", "(@1*cos(@0) - @2*sin(@0))/(sin(@0)+cos(@0))",
                         ROOT.RooArgList(theta, m_x, m_y))

s_u = ROOT.RooFormulaVar("#sigma_{u}", "#sigma_{u}",
                         "sqrt((sin(@0)*@2)**2+(cos(@0)*@3)**2+2*sin(@0)*cos(@0)*@1*@2*@3)/(sin(@0)+cos(@0))",
                         ROOT.RooArgList(theta, rho, s_x, s_y))
s_v = ROOT.RooFormulaVar("#sigma_{v}", "#sigma_{v}",
                         "sqrt((cos(@0)*@2)**2+(sin(@0)*@3)**2-2*sin(@0)*cos(@0)*@1*@2*@3)/(sin(@0)+cos(@0))",
                         ROOT.RooArgList(theta, rho, s_x, s_y))

E_u = ROOT.RooFormulaVar("E_u", "u Energy", "(@1*sin(@0) + @2*cos(@0))/(sin(@0)+cos(@0))",
                         ROOT.RooArgList(theta, E_x, E_y))
E_v = ROOT.RooFormulaVar("E_v", "v Energy", "(@1*cos(@0) - @2*sin(@0))/(sin(@0)+cos(@0))",
                         ROOT.RooArgList(theta, E_x, E_y))
E_u_plot = ROOT.RooRealVar("E_u", "Rotated Energy", E_lo, E_hi, "E/E_{0}")

g = ROOT.TGraph()
s_min = m_u.getVal()
old_th = theta.getVal()
th_min = 0
for test_th in numpy.linspace(0.75*theta.getVal(), 1.25*theta.getVal(), 1000, False):
    n = g.GetN()
    g.Set(n+1)
    theta.setVal(test_th)
    s = s_u.getVal()/E_u.getVal()
    if s < s_min:
        s_min = s
        th_min = test_th
    g.SetPoint(n, test_th*180/pi, s)
theta.setVal(old_th)

print("Min sigma/E at  theta = %0.2f"%(th_min*180/pi))

c_g = ROOT.TCanvas("c_g")
g.Draw("AP")
g.GetXaxis().SetTitle("#theta (deg)")
g.GetYaxis().SetTitle("#sigma_{rot}/E")
    

f_x = ROOT.RooGaussian("f_x", "f_x", E_x, m_x, s_x)
f_y= ROOT.RooGaussian("f_y", "f_y", E_y, m_y, s_y)

f_u = ROOT.RooGaussian("f_u", "f_u", E_u, m_u, s_u)
f_v = ROOT.RooGaussian("f_v", "f_v", E_v, m_v, s_v)
f_2d = ROOT.RooProdPdf("f_2d", "f_2d", f_u, f_v)
f_u_plot = ROOT.RooGaussian("f_u_plot", "f_u_plot", E_u_plot, m_u, s_u)


h_2d = E_x.createHistogram("h_2d", E_y)
f_2d.fillHistogram(h_2d, ROOT.RooArgList(E_x, E_y))

t_x = ROOT.TPaveText(0.70, 0.90, 0.98, 0.98, "ndc")
t_x.SetBorderSize(0)
t_x.SetFillColor(ROOT.kWhite)
t_x.AddText("#sigma/E = %0.1f%%"%(100*s_x.getVal()/m_x.getVal()))

t_y = ROOT.TPaveText(0.70, 0.90, 0.98, 0.98, "ndc")
t_y.SetBorderSize(0)
t_y.SetFillColor(ROOT.kWhite)
t_y.AddText("#sigma/E = %0.1f%%"%(100*s_y.getVal()/m_y.getVal()))

t_u = ROOT.TPaveText(0.70, 0.90, 0.98, 0.98, "ndc")
t_u.SetBorderSize(0)
t_u.SetFillColor(ROOT.kWhite)
t_u.AddText("#sigma/E = %0.1f%%"%(100*s_u.getVal()/m_u.getVal()))

ax_x_lo = E_lo
ax_y_lo = E_lo
ax_w_lo = E_lo

ax_x_hi = min(E_hi, ax_x_lo + (E_hi-ax_y_lo)*tan(theta.getVal()))
ax_y_hi = min(E_hi, ax_y_lo + (E_hi-ax_x_lo)/tan(theta.getVal()))
ax_w_hi = E_lo + ((ax_y_hi - ax_y_lo)*cos(theta.getVal())+(ax_x_hi-ax_x_lo)*sin(theta.getVal()))/(sin(theta.getVal())+cos(theta.getVal()))
#print("(%0.2f %0.2f) (%0.2f %0.2f) %0.2f %0.2f"%(ax_x_lo, ax_y_lo, ax_x_hi, ax_y_hi, ax_w_lo, ax_w_hi))
ax = ROOT.TGaxis(ax_x_lo, ax_y_lo, ax_x_hi, ax_y_hi, ax_w_lo, ax_w_hi, 410, "-")
ax.CenterTitle(True)
ax.SetTitleOffset(1.5)
ax.SetLabelOffset(0.05)
ax.SetTitleColor(ROOT.kGray+2)
ax.SetLabelColor(ROOT.kGray+2)
ax.SetLineColor(ROOT.kGray+2)
ax.SetTitle("Rotated Energy (E/E_{0})")

c_x = ROOT.TCanvas("c_x", "c_x", 600, 600)
c_x.SetRightMargin(0.01)
c_x.SetTopMargin(0.01)
frame_x = E_x.frame()
f_x.plotOn(frame_x)
frame_x.Draw()
t_x.Draw()
frame_x.GetYaxis().SetTitle("dN/dE")
frame_x.GetYaxis().SetTitleOffset(1.25)
frame_x.SetTitle("")
c_x.Print("toy_anticorrelation_scint.pdf")

c_y = ROOT.TCanvas("c_y", "c_y", 600, 600)
c_y.SetRightMargin(0.01)
c_y.SetTopMargin(0.01)
frame_y = E_y.frame()
f_y.plotOn(frame_y)
frame_y.Draw()
t_y.Draw()
frame_y.GetYaxis().SetTitle("dN/dE")
frame_y.GetYaxis().SetTitleOffset(1.25)
frame_y.SetTitle("")
c_y.Print("toy_anticorrelation_ioniz.pdf")

c_u = ROOT.TCanvas("c_u", "c_u", 600, 600)
c_u.SetRightMargin(0.01)
c_u.SetTopMargin(0.01)
frame_u = E_u_plot.frame()
f_u_plot.plotOn(frame_u)
frame_u.Draw()
t_u.Draw()
frame_u.GetYaxis().SetTitle("dN/dE")
frame_u.GetYaxis().SetTitleOffset(1.25)
frame_u.SetTitle("")
c_u.Print("toy_anticorrelation_rot.pdf")

c_2d = ROOT.TCanvas("c_2d", "c_2d", 600, 600)
c_2d.SetRightMargin(0.01)
c_2d.SetTopMargin(0.01)
h_2d.Draw("CONT0")
ax.Draw()
h_2d.SetTitle("")
c_2d.Print("toy_anticorrelation_2d.pdf")

c4 = ROOT.TCanvas("c4", "c4", 600, 600)
c4.Divide(2, 2)

c4.cd(1)
frame_y.Draw()
t_y.Draw()

c4.cd(2)
h_2d.Draw("CONT0")
ax.Draw()

c4.cd(3)
frame_u.Draw()
t_u.Draw()

c4.cd(4)
frame_x.Draw()
t_x.Draw()



## c4.cd(2)
## f_2d.Draw("CONT")
## #ROOT.gPad.SetLogz()
## c4.cd(3)
## f_1d.Draw()
## c4.cd(4)
## f_x.Draw()

## c2 = ROOT.TCanvas("c2", "c2", 1280, 480)
## c2.Divide(2, 1)

## c2.cd(1)
## f_2d.Draw("CONT")
## #ROOT.gPad.SetLogz()

## c2.cd(2)
## f_1d.Draw()
## f_x.Draw("same")
## f_y.Draw("same")

dummy = raw_input("Press enter...")
