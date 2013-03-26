#!/bin/env python

import ROOT
import math
import array

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
theta = ROOT.RooRealVar("#theta", "#theta", math.atan2(s_y.getVal(), s_x.getVal()))
print("theta = %0.2f"%(theta.getVal()*180/math.pi))

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
h_2d.SetTitle("")
c_2d.Print("toy_anticorrelation_2d.pdf")

c4 = ROOT.TCanvas("c4", "c4", 600, 600)
c4.Divide(2, 2)

c4.cd(1)
frame_y.Draw()
t_y.Draw()

c4.cd(2)
h_2d.Draw("CONT0")

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