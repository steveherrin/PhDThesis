#!/bin/env python

import array
from math import pi, sin, cos, asin, acos, fabs, sqrt
import ROOT
ROOT.gSystem.Load("libEXOUtilities")
ROOT.gStyle.SetOptStat(0)
ROOT.TGaxis.SetMaxDigits(3)

muon = ROOT.TChain("muon")
muon.Add("simulations/muonsim_*.root")

use_symmetry = False

max_theta = pi/2

def SetBins(h):
    # Set the binning on a TH2Poly in a way that's useful for zenith
    # angle stuff. The bin with the zenith is not divided in phi. The
    # next bin in theta is divided into n_p_top bins. Subsequent theta
    # bins are divided into phi bins with a similar solid angle.

    # n_t = 10, n_p = 18, n_p_top = 4 shows off the wire plane holes well
    
    n_t = 90 # n theta bins
    n_p_top = 12 # 4-fold symmetry to start with
    
    if use_symmetry:
        sa_top = pi/2*(cos(max_theta/n_t)-cos(2*max_theta/n_t))
    else:
        sa_top = 2*pi*(cos(max_theta/n_t)-cos(2*max_theta/n_t))
    desired_sa = sa_top/n_p_top
        
    for i in xrange(n_t):
        #if i == 0:
        #    t_lo = 0
        #    t_hi = pi/180
        #else:
        t_lo = (max_theta)*(i)/(n_t)
        t_hi = (max_theta)*(i+1)/(n_t)

        if use_symmetry:
            sa = pi/2*(cos(t_lo)-cos(t_hi))
        else:
            sa = 2*pi*(cos(t_lo)-cos(t_hi))

        if (i == 0):
            n_p = 1
        else:
            # in case we want to find even-numbered divisions
            d = 4
            n_p = d
            while (sa/n_p > desired_sa):
                n_p += d
            # if we divided too finely
            if (n_p > d
                and fabs(sa/(n_p-d)-desired_sa) < fabs(sa/(n_p)-desired_sa)):
                n_p -= d
        n_p = 120
            
        
        for j in xrange(n_p):
            if use_symmetry:
                p_lo = j*pi/(2*n_p)
                p_hi = (j+1)*pi/(2*n_p)
            else:
                p_lo = 2*j*pi/(n_p)-pi
                p_hi = 2*(j+1)*pi/(n_p)-pi
            
            
            h.AddBin(p_lo, t_lo, p_hi, t_hi)

def Add(h1, h2, c1=1, c2=1):
    # Add h1 and h2 and return a new histogram with the result
    h = h1.Clone(h1.GetName()+"_plus_"+h2.GetName())
    n = h1.fN
    if (n != h2.fN):
        raise ValueError("Histograms don't have same binning.")
    for i in xrange(0, n):
        h.SetBinContent(i, c1*h1.GetBinContent(i)+c2*h2.GetBinContent(i))
    return h
    
def Divide(h1, h2, c1=1):
    # Divide h1 by h2 and return a new histogram with the result
    h = h1.Clone(h1.GetName()+"_div_"+h2.GetName())
    n = h1.fN
    if (n != h2.fN):
        raise ValueError("Histograms don't have same binning.")
    for i in xrange(0, n):
        try:
            h.SetBinContent(i, c1*float(h1.GetBinContent(i))/h2.GetBinContent(i))
        except ZeroDivisionError:
            h.SetBinContent(i, 0)
    return h

def SqrtAndDivide(h1, h2, c1=1):
    # Take sqrt of h1 and divide by h2 and return a new histogram with the result
    h = h1.Clone(h1.GetName()+"_div_"+h2.GetName())
    n = h1.fN
    if (n != h2.fN):
        raise ValueError("Histograms don't have same binning.")
    for i in xrange(0, n):
        try:
            h.SetBinContent(i, c1*sqrt(h1.GetBinContent(i))/h2.GetBinContent(i))
        except ZeroDivisionError:
            h.SetBinContent(i, 0)
    return h
            
def ProjectedArea(theta, phi):
    r = 22.875-0.137
    h = 2*20.44065
    if use_symmetry:
        n = 4
    else:
        n = 1
    return n*((pi*r**2*fabs(cos(phi))*sin(theta) +
               2*r*h*sqrt(1-sin(theta)**2*cos(phi)**2)))

def AddToHist(h, theta, phi, to_add):
    # Add to_add to whatever is already in bin theta, phi
    bin_n = h.FindBin(phi, theta)
    old = h.GetBinContent(bin_n)
    new = old + to_add
    h.SetBinContent(bin_n, to_add)

def KeepMax(h, theta, phi, value):
    # replace what's in bin theta, phi with value if value is greater
    # than whatever's already there

    bin_n = h.FindBin(phi, theta)
    old = h.GetBinContent(bin_n)
    if value > old:
        h.SetBinContent(bin_n, value)

def MakeZeroNeg(h1, h2):
    # Make what's in h1 negative if h2 is zero
     n = h1.fN
     if (n != h2.fN):
        raise ValueError("Histograms don't have same binning.")
     for i in xrange(0, n):
         if h2.GetBinContent(i) == 0:
             h1.SetBinContent(i, -1)


#t_bins = array.array('d',[i*pi/32 for i in xrange(0,17)])
#p_bins = array.array('d',[pi*(0.125*i-1) for i in xrange(0,17)])

if use_symmetry:
    p_lo = 0
    p_hi = pi/2
else:
    p_lo = -pi
    p_hi = pi

#n_p = 120
#n_t = 90
n_p = 30
n_t = 20

h_sim = ROOT.TH2D("h_sim","Number Simulated",
             n_p, p_lo, p_hi, n_t, 0, max_theta)

h_corr = ROOT.TH2D("h_corr","Number Correctly Reconstructed",
                      n_p, p_lo, p_hi, n_t, 0, max_theta)

h_incorr = ROOT.TH2D("h_incorr","Number Incorrectly Reconstructed",
                        n_p, p_lo, p_hi, n_t, 0, max_theta)

h_acc = ROOT.TH2D("h_acc","Acceptance",
                     n_p, p_lo, p_hi, n_t, 0, max_theta)

h_mask = ROOT.TH2D("h_mask","Mask",
                     n_p, p_lo, p_hi, n_t, 0, max_theta)

h_terrc = ROOT.TH2D("h_terrc","Theta Error Sum",
                      n_p, p_lo, p_hi, n_t, 0, max_theta)

h_perrc = ROOT.TH2D("h_perrc","Phi Error Sum",
                      n_p, p_lo, p_hi, n_t, 0, max_theta)

h_tperrc = ROOT.TH2D("h_tperrc","Combined Error Sum",
                      n_p, p_lo, p_hi, n_t, 0, max_theta)

h_terrm = ROOT.TH2D("h_terrm","Theta Error Max",
                       n_p, p_lo, p_hi, n_t, 0, max_theta)

h_perrm = ROOT.TH2D("h_perrm","Phi Error Max",
                      n_p, p_lo, p_hi, n_t, 0, max_theta)

h_tperrm = ROOT.TH2D("h_tperrm","Combined Error Max",
                      n_p, p_lo, p_hi, n_t, 0, max_theta)

h_lterr = ROOT.TH2D("h_lterr","Number Large Theta Error",
                       n_p, p_lo, p_hi, n_t, 0, max_theta)
h_lperr = ROOT.TH2D("h_lperr","Number Large Phi Error",
                       n_p, p_lo, p_hi, n_t, 0, max_theta)
h_ltperr = ROOT.TH2D("h_ltperr","Number Large Errors",
                       n_p, p_lo, p_hi, n_t, 0, max_theta)



for h in [h_sim, h_corr, h_incorr, h_acc, h_terrc, h_perrc, h_tperrc, h_terrm, h_perrm, h_tperrm, h_mask, h_lterr, h_lperr, h_ltperr]:
    #SetBins(h)
    h.GetXaxis().SetTitle("#phi")
    h.GetYaxis().SetTitle("#theta")

md = ROOT.EXOMuonData()
muon.SetBranchAddress("MuonBranch", md)

N_reconstructed = 0
N_ismuon = 0
N_othercheck = 0

n = muon.GetEntries()
for i in xrange(n):
    if (i %(n//10) == 0):
        print("%i of %i."%(i, n))
    muon.GetEntry(i)

    if md.ContainsMuon():
        N_ismuon += 1

    weight = 0
    for j in xrange(md.GetNumMuonTracks()):
        track = md.GetMuonTrack(j)

        if track.fTheta >= 0:
            weight += 1

    if weight > 0:
        N_othercheck += 1
    elif md.ContainsMuon():
        print("WTF? Contains muon but weight = %f"%(weight))
    

    if weight <= 0:
        continue

    mc_theta = md.fMonteCarloTheta

    if use_symmetry:
        if (md.fMonteCarloPhi < 0):
            mc_phi = fabs(-pi/2 - md.fMonteCarloPhi)
        else:
            mc_phi = fabs(pi/2 - md.fMonteCarloPhi)
    else:
        mc_phi = md.fMonteCarloPhi

    h_sim.Fill(mc_phi, mc_theta)

    true_bin_n = h_sim.FindBin(mc_phi, md.fMonteCarloTheta)

    for j in xrange(md.GetNumMuonTracks()):
        track = md.GetMuonTrack(j)

        theta = track.fTheta

        if use_symmetry:
            if (track.fPhi < 0):
                phi = fabs(-pi/2 - track.fPhi)
            else:
                phi = fabs(-pi/2 - track.fPhi)
        else:
            phi = track.fPhi

        if (phi > -900 and theta > -900):
            N_reconstructed += 1.0/weight
            

        recon_bin = h_sim.FindBin(phi, theta)

        if (recon_bin == true_bin_n):
            h_corr.Fill(phi, theta, 1.0/weight)
        else:
            h_incorr.Fill(phi, theta, 1.0/weight)
            #h_incorr.Fill(md.fMonteCarloPhi, md.fMonteCarloTheta,
            #              1./md.GetNumMuonTracks())

            phi_err = acos(cos(phi - mc_phi))
            #theta_err = sin(mc_theta)*acos(cos(theta - mc_theta))
            theta_err = acos(cos(theta - mc_theta))
            tp_err = 2*asin(sqrt(sin((theta-mc_theta)/2)**2
                                 +cos(theta)*cos(mc_theta)
                                 *sin((phi-mc_phi)/2)**2))/(2*pi)

            if (theta_err > 5*pi/180):
                h_lterr.Fill(phi, theta, 1.0/weight)
            if (phi_err > 5*pi/180):
                h_lperr.Fill(phi, theta, 1.0/weight)
            if (tp_err > 5*pi/180):
                h_ltperr.Fill(phi, theta, 1.0/weight)
        
            #AddToHist(h_perrc, theta, phi, phi_err**2)
            #AddToHist(h_terrc, theta, phi, theta_err**2)
            #AddToHist(h_tperrc, theta, phi, tp_err**2)
 
            #KeepMax(h_perrm, theta, phi, phi_err)
            #KeepMax(h_terrm, theta, phi, theta_err)
            #KeepMax(h_tperrm, theta, phi, tp_err)

#good_phi = [0, 0.15, 0.9, 1.15, 1.2, 0.8, 0]
# below is good 2013-03-01
# good_phi = [0, 0.15, 0.9, 1.2, 1.2, 0.1]
good_phi = [0, 0.15, 0.9, 1.2, 1.2, 0.8]
#good_phi.extend([-x for x in reversed(good_phi[:-1])])
good_phi.extend([-x for x in reversed(good_phi)])
good_phi.append(good_phi[0])
#good_theta = [0.4, 0.45, 0.45, 0.75, 0.95, 1, 1.3]
# below is good 2013-03-1
#good_theta = [0.4, 0.45, 0.45, 0.75, 0.95, 1.3]
good_theta = [0.4, 0.45, 0.45, 0.75, 0.95, pi/3]
#good_theta.extend(reversed(good_theta[:-1]))
good_theta.extend(reversed(good_theta))
good_theta.append(good_theta[0])

g_good_region_pos = ROOT.TGraph(len(good_phi))
g_good_region_neg = ROOT.TGraph(len(good_phi))
for i in xrange(len(good_phi)):
    g_good_region_pos.SetPoint(i, pi-good_phi[i], good_theta[i])
    g_good_region_neg.SetPoint(i, good_phi[i]-pi, good_theta[i])
g_good_region_pos.SetLineColor(ROOT.kBlack)
g_good_region_neg.SetLineColor(ROOT.kBlack)

good_phi.extend([-x for x in reversed(good_phi[:-1])])
good_theta.extend(reversed(good_theta[:-1]))

good_phi_a = array.array('d', good_phi)
good_theta_a = array.array('d', good_theta)
g_good_region = ROOT.TGraph(len(good_phi_a), good_phi_a, good_theta_a)
g_good_region.SetLineColor(ROOT.kBlack)
        

h_cplusi = Add(h_corr, h_incorr)
h_cplusi.SetNameTitle("h_cplusi","Number Reconstructed")
            
h_eff = Divide(h_corr, h_sim)
h_eff.SetNameTitle("h_eff","Efficiency (correct only)")

h_mis = Divide(h_incorr, h_sim)
h_mis.SetNameTitle("h_mis","Misreconstructed Rate")

h_ratio = Divide(h_corr, h_incorr)
h_ratio.SetNameTitle("h_ratio","(Correctly Reconstructed Rate)/(Misreconstructed Rate)")

h_eff2 = Divide(h_cplusi, h_sim)
h_eff2.SetNameTitle("h_eff2","Efficiency (all reconstructed)")

h_terrrate = Divide(h_lterr, h_cplusi)
print("h_terrrate min %f"%(h_terrrate.GetMinimum()))
MakeZeroNeg(h_terrrate, h_cplusi)
h_terrrate.SetMinimum(-1e-2)
h_terrrate.SetNameTitle("h_terrrate","Rate of large #theta errors")

h_perrrate = Divide(h_lperr, h_cplusi)
print("h_perrrate min %f"%(h_perrrate.GetMinimum()))
MakeZeroNeg(h_perrrate, h_cplusi)
h_perrrate.SetMinimum(-1e-2)
h_perrrate.SetNameTitle("h_terrrate","Rate of large #phi errors")

h_tperrrate = Divide(h_ltperr, h_cplusi)
print("h_tperrrate min %f"%(h_tperrrate.GetMinimum()))
MakeZeroNeg(h_tperrrate, h_cplusi)
h_tperrrate.SetMinimum(-1e-2)
h_tperrrate.SetNameTitle("h_tperrrate","Rate of large angular errors")

for i in xrange(h_acc.fN):
    x = ROOT.Long(0)
    y = ROOT.Long(0)
    z = ROOT.Long(0)
    h_acc.GetBinXYZ(i, x, y, z)
    theta = (h_acc.GetYaxis().GetBinLowEdge(y) + h_acc.GetYaxis().GetBinUpEdge(y))/2
    phi = (h_acc.GetXaxis().GetBinLowEdge(x) + h_acc.GetXaxis().GetBinUpEdge(x))/2

    eff = h_eff2.GetBinContent(i)
    h_acc.SetBinContent(i, eff*ProjectedArea(theta, phi))
    

#h_perr = SqrtAndDivide(h_perrc, h_incorr)
#h_perr.SetNameTitle("h_perr","RMS Phi Error")

#h_terr = SqrtAndDivide(h_terrc, h_incorr)
#h_terr.SetNameTitle("h_terr","RMS Theta Error")

#h_tperr = SqrtAndDivide(h_tperrc, h_incorr)
#h_tperr.SetNameTitle("h_tperr","RMS Combined Error")

h_err = Divide(h_incorr, h_cplusi)
#h_err = Divide(h_incorr, h_sim)
h_err.SetNameTitle("h_err","Error Rate")

c_eff = ROOT.TCanvas("c_eff")
h_eff.Draw("COLZ")
g_good_region.Draw("L")
g_good_region_pos.Draw("L")
g_good_region_neg.Draw("L")
c_mis = ROOT.TCanvas("c_mis")
h_mis.Draw("COLZ")
c_ratio = ROOT.TCanvas("c_ratio")
h_ratio.Draw("COLZ")
#c_mask = ROOT.TCanvas("c_mask")
#h_mask.Draw("COLZ")
c_eff2 = ROOT.TCanvas("c_eff2")
h_eff2.Draw("COLZ")
g_good_region.Draw("L")
g_good_region_pos.Draw("L")
g_good_region_neg.Draw("L")
c_terrrate = ROOT.TCanvas("c_terrrate")
h_terrrate.Draw("COLZ")
g_good_region.Draw("L")
g_good_region_pos.Draw("L")
g_good_region_neg.Draw("L")
c_perrrate = ROOT.TCanvas("c_perrrate")
h_perrrate.Draw("COLZ")
g_good_region.Draw("L")
g_good_region_pos.Draw("L")
g_good_region_neg.Draw("L")
c_tperrrate = ROOT.TCanvas("c_tperrrate")
h_tperrrate.Draw("COLZ")
g_good_region.Draw("L")
g_good_region_pos.Draw("L")
g_good_region_neg.Draw("L")
#c_perr = ROOT.TCanvas("c_perr")
#h_perr.Draw("COLZ")
#c_terr = ROOT.TCanvas("c_terr")
#h_terr.Draw("COLZ")
#c_tperr = ROOT.TCanvas("c_tperr")
#h_tperr.Draw("COLZ")
#c_perrm = ROOT.TCanvas("c_perrm")
#h_perrm.Draw("COLZ")
#c_terrm = ROOT.TCanvas("c_terrm")
#h_terrm.Draw("COLZ")
#c_tperrm = ROOT.TCanvas("c_tperrm")
#h_tperrm.Draw("COLZ")
c_err = ROOT.TCanvas("c_err")
h_err.Draw("COLZ")
c_acc = ROOT.TCanvas("c_acc")
h_acc.Draw("COLZ")

print("Out of %i simulated, tagged %f (= %0.2f) %%"%(n, N_reconstructed,
                                                     100.0*N_reconstructed/n))
print("Out of %i simulated, tagged %f (= %0.2f) %%"%(n, N_ismuon,
                                                     100.0*N_ismuon/n))
print("Out of %i simulated, tagged %f (= %0.2f) %%"%(n, N_othercheck,
                                                     100.0*N_othercheck/n))

dummy = raw_input("Press Enter...")
