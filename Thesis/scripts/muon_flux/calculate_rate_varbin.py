#!/bin/env python

import array
from math import pi, sin, cos, tan, acos, fabs, sqrt, pow, log, exp
import scipy.integrate
import glob
import ROOT
ROOT.gSystem.Load("libEXOUtilities")
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(111)
ROOT.TGaxis.SetMaxDigits(3)

use_symmetry = False
max_theta = 80*pi/180

def flux(theta, h):

    log_intensity = (-(8e-4*h)/cos(theta)+(8e-4*h))
    
    intensity = exp(log_intensity)
    intensity *= pow(cos(theta), 1.53)
    return intensity

def flux_sin(theta, h):
    intensity = flux(theta, h)*sin(theta)
    return intensity


def vertical_flux_sin(theta, h):
    intensity = flux(theta, h)*cos(theta)*sin(theta)
    return intensity

def GetFluxRatio(theta_lo, theta_hi, phi_lo, phi_hi):
    depth = 1585

    ratio = ((phi_hi-phi_lo)*scipy.integrate.quad(flux_sin,
                                                  theta_lo, theta_hi,
                                                  (depth,))[0]
             /(2*pi*scipy.integrate.quad(flux_sin, 0, pi/2, (depth))[0]))
    return ratio

def GetVFluxRatio(theta_lo, theta_hi, phi_lo, phi_hi):
    depth = 1585

    ratio = ((phi_hi-phi_lo)*scipy.integrate.quad(flux_sin,
                                                  theta_lo, theta_hi,
                                                  (depth,))[0]
             /(2*pi*scipy.integrate.quad(vertical_flux_sin, 0, pi/2, (depth))[0]))
    return ratio

def GetGoldenTime():
    dT = 0
    run_list = open("golden_runs_20130213.txt")
    for run in run_list:
        try:
            files = glob.glob("golden_runs/muon_%s*.root"%(run.rstrip()))
            files.sort()

            f = ROOT.TFile(files.pop())
            record_list = f.Get("tree").GetUserInfo().At(1)

            begin = record_list.GetNextRecord('EXOBeginRecord')()
            begin_t = begin.GetTimestamp()
            end = record_list.GetNextRecord('EXOEndRecord')()
            end_t = end.GetTimestamp()

        except (AttributeError, IndexError, ReferenceError):
            print("Problem with run %s."%(run.rstrip()))
            continue

        dT += end_t - begin_t

    return dT/(1e9) # so it's in seconds
        
def SetBins(h):
    # Set the binning on a TH2Poly in a way that's useful for zenith
    # angle stuff. The bin with the zenith is not divided in phi. The
    # next bin in theta is divided into n_p_top bins. Subsequent theta
    # bins are divided into phi bins with a similar solid angle.

    # n_t = 10, n_p = 18, n_p_top = 4 shows off the wire plane holes well
    
    n_t = 80//1 # n theta bins
    n_p_top = 8 # 4-fold symmetry to start with
    
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
        n_p = 360//3
            
        
        for j in xrange(n_p):
            if use_symmetry:
                p_lo = j*pi/(2*n_p)
                p_hi = (j+1)*pi/(2*n_p)
            else:
                p_lo = 2*j*pi/(n_p)-pi
                p_hi = 2*(j+1)*pi/(n_p)-pi
            
            
            h.AddBin(p_lo, t_lo, p_hi, t_hi)    

def Collapse(h):
    h_1d = ROOT.TH1D(h.GetName()+"_bins", h.GetTitle() + " Bin Distribution",
                     100, 0, 1e-6)
    for bin in h.GetBins():
        if (bin.GetContent() > 0):
            h_1d.Fill(bin.GetContent())
    h_1d.GetXaxis().SetTitle("Vertical Flux (cm^{-2} sr^{-1} s^{-1})")
    return h_1d

def Zero(h):
    for i in xrange(h.GetNumberOfBins()):
        h.SetBinContent(i+1, 0)

def GetFluxAndError(N, N_rec, N_sim, denominator):
    if N == 0:
        return (0, 0)
    try:
        eff = N_rec/N_sim
    except ZeroDivisionError:
        print("eff = ; N_sim = %f"%(N_sim))
    eff_err = sqrt(N_rec)/N_sim
    try:
        N_err = sqrt(N)/eff
    except ZeroDivisionError:
        print("N_err = ; eff = %f"%(eff))
    try:
        flux = N/(eff*denominator)
    except ZeroDivisionError:
        print("flux = ; denominator = %f"%(denominator))
    try:
        flux_err = flux*sqrt((N_err/N)**2 + (eff_err/eff)**2)
    except ZeroDivisionError:
        print("flux_err = ; N = %f, eff = %f"%(N, eff))
    return (flux, flux_err)

def CollapseTheta(h):
    edges = [0, max_theta]
    
    for bin in h.GetBins():
        edge = bin.GetYMin()
        if edge not in edges:
            edges.append(edge)
    edges.sort()
    nums = [1 for i in xrange(len(edges)-1)]
    edges_array = array.array('d', edges)
    h_1d = ROOT.TH1D(h.GetName()+"_py", h.GetTitle() + " #theta projection",
                     len(edges)-1, edges_array)
    for bin in h.GetBins():
        theta = bin.GetYMin()
        n = bin.GetContent()
        h_1d.Fill(theta, n)
        nums[h_1d.FindBin(theta)-1] += 1
    for i, num in enumerate(nums):
        x = h_1d.GetBinContent(i+1)
        h_1d.SetBinContent(i+1, x/num)
    h_1d.GetXaxis().SetTitle("#theta")
    return h_1d
    

def FillHistFromTreeMC(h, chain):
    md = ROOT.EXOMuonData()
    chain.SetBranchAddress("MuonBranch", md)

    n = chain.GetEntries()
    print("Processing %i events."%(n))
    for i in xrange(n):
        if (i > 0 and i % (n//5) == 0):
            print("  %i of %i..."%(i,n))

        chain.GetEntry(i)

        mc_theta = md.fMonteCarloTheta

        if use_symmetry:
            if (md.fMonteCarloPhi < 0):
                mc_phi = fabs(-pi/2 - md.fMonteCarloPhi)
            else:
                mc_phi = fabs(pi/2 - md.fMonteCarloPhi)
        else:
            mc_phi = md.fMonteCarloPhi

        if (mc_theta < 0 or fabs(mc_phi) > 2*pi):
            continue

        h.Fill(mc_phi, mc_theta)
    

def FillHistFromTree(h, chain):
    md = ROOT.EXOMuonData()
    chain.SetBranchAddress("MuonBranch", md)

    n = chain.GetEntries()
    print("Processing %i events."%(n))
    for i in xrange(n):
        if (i > 0 and i % (n//5) == 0):
            print("  %i of %i..."%(i,n))

        chain.GetEntry(i)

        weight = 0
        for j in xrange(md.GetNumMuonTracks()):
            track = md.GetMuonTrack(j)

            if track.fTheta >= 0:
                weight += 1
 
        for j in xrange(md.GetNumMuonTracks()):
            track = md.GetMuonTrack(j)

            #if (track.fNumUHitsOnTrack < 4 or track.fNumVHitsOnTrack < 4 or
            #    track.fNumUHitsOnTrack + track.fNumUHitsOnTrack < 10):
            #    continue

            theta = track.fTheta

            if use_symmetry:
                if (track.fPhi < 0):
                    phi = fabs(-pi/2 - track.fPhi)
                else:
                    phi = fabs(-pi/2 - track.fPhi)
            else:
                phi = track.fPhi

            if (theta < 0 or fabs(phi) > pi):
                continue
    
            h.Fill(phi, theta, 1./weight)


    

def Add(h1, h2, c1=1, c2=1):
    # Add h1 and h2 and return a new histogram with the result
    h = h1.Clone(h1.GetName()+"_plus_"+h2.GetName())
    n = h1.GetNumberOfBins()
    if (n != h2.GetNumberOfBins()):
        raise ValueError("Histograms don't have same binning.")
    for i in xrange(1, n+1):
        h.SetBinContent(i, c1*h1.GetBinContent(i)+c2*h2.GetBinContent(i))
    return h
    
def Divide(h1, h2, c1=1):
    # Divide h1 by h2 and return a new histogram with the result
    h = h1.Clone(h1.GetName()+"_div_"+h2.GetName())
    n = h1.GetNumberOfBins()
    if (n != h2.GetNumberOfBins()):
        raise ValueError("Histograms don't have same binning.")
    for i in xrange(1, n+1):
        try:
            h.SetBinContent(i, c1*float(h1.GetBinContent(i))/h2.GetBinContent(i))
        except ZeroDivisionError:
            h.SetBinContent(i, 0)
    return h

def SqrtAndDivide(h1, h2, c1=1):
    # Take sqrt of h1 and divide by h2 and return a new histogram with the result
    h = h1.Clone(h1.GetName()+"_div_"+h2.GetName())
    n = h1.GetNumberOfBins()
    if (n != h2.GetNumberOfBins()):
        raise ValueError("Histograms don't have same binning.")
    for i in xrange(1, n+1):
        try:
            h.SetBinContent(i, c1*sqrt(h1.GetBinContent(i))/h2.GetBinContent(i))
        except ZeroDivisionError:
            h.SetBinContent(i, 0)
    return h
            
def ProjectedArea(theta, phi):
    # projected area for a particle incident from theta, phi
    # in cm^2
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

#t_bins = array.array('d',[i*pi/32 for i in xrange(0,17)])
#p_bins = array.array('d',[pi*(0.125*i-1) for i in xrange(0,17)])

if use_symmetry:
    p_lo = 0
    p_hi = pi/2
else:
    p_lo = -pi
    p_hi = pi

h_sim = ROOT.TH2Poly("h_sim","Number Simulated",
             p_lo, p_hi, 0, max_theta)

h_rec = ROOT.TH2Poly("h_rec","Number Reconstructed",
                     p_lo, p_hi, 0, max_theta)

h_obs = ROOT.TH2Poly("h_obs","Number Observed",
                     p_lo, p_hi, 0, max_theta)


for h in [h_sim, h_rec, h_obs]:
    SetBins(h)
    h.GetXaxis().SetTitle("#phi")
    h.GetYaxis().SetTitle("#theta")

mc_muon = ROOT.TChain("muon")
mc_muon.Add("simulations/muonsim_*.root")
FillHistFromTreeMC(h_sim, mc_muon)
FillHistFromTree(h_rec, mc_muon)

muon = ROOT.TChain("muon")
muon.Add("golden_runs/muon_*.root")
FillHistFromTree(h_obs, muon)
            
h_eff = Divide(h_rec, h_sim)
h_eff.SetNameTitle("h_eff","Efficiency")

h_num = Divide(h_obs, h_eff)
h_num.SetNameTitle("h_num","Observed Number (corrected for efficiency)")

h_flux = h_num.Clone()
Zero(h_flux)
h_flux.SetNameTitle("h_flux","Muon Flux")

h_vflux = h_num.Clone()
Zero(h_vflux)
h_vflux.SetNameTitle("h_vflux","Vertical Muon Flux")

h_check = h_num.Clone()
Zero(h_check)
h_check.SetNameTitle("h_check","Region Check")

#g_vflux = ROOT.TGraphErrors()
#g_flux = ROOT.TGraphErrors()
#g_flux_v_sec = ROOT.TGraphErrors()
#g_vflux_scaled = ROOT.TGraphErrors()
#g_flux_scaled = ROOT.TGraphErrors()

h_I = ROOT.TH1D("h_I", "I(h, 0)", 100, 0, 1e-6)


dT = GetGoldenTime()
N_sim = 0
N_rec = 0
N_total = 0
sum_flux = 0
sum_vflux = 0
sum_sq_flux_err = 0
sum_sq_vflux_err = 0
sum_sa = 0
#first_bin_theta = 0

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

good_phi_a = array.array('d', good_phi)
good_theta_a = array.array('d', good_theta)
g_good_region = ROOT.TGraph(len(good_phi_a), good_phi_a, good_theta_a)

for i, bin in enumerate(h_flux.GetBins()):
    if bin.GetYMin() == 0:
        theta = 0
        #first_bin_theta = bin.GetYMax()
        #h_flux.SetBinContent(i+1, 0)
        #h_vflux.SetBinContent(i+1, 0)
        #continue
    else:
        #theta = (bin.GetYMin()+bin.GetYMax())/2
        #theta = bin.GetYMax()
        theta = (scipy.integrate.quad(lambda x: x*sin(x), bin.GetYMin(),
                                      bin.GetYMax())[0]/
                 scipy.integrate.quad(lambda x: sin(x), bin.GetYMin(),
                                      bin.GetYMax())[0])
    
    phi = (bin.GetXMin()+bin.GetXMax())/2

    if (phi < -pi/2):
        phi_sym = phi + pi
    elif (phi > pi/2):
        phi_sym = phi - pi
    else:
        phi_sym = phi

    #h_check.SetBinContent(i+1, g_good_region.IsInside(phi_sym, theta))
    
        

    if (g_good_region.IsInside(phi_sym, theta) == 1):

        # solid angle of the bin
        d_phi = bin.GetXMax() - bin.GetXMin()
        d_theta = bin.GetYMax() - bin.GetYMin()
        sa = d_phi*(cos(bin.GetYMin())-cos(bin.GetYMax()))
        sum_sa += sa
    
        N = h_num.GetBinContent(i+1)
        h_check.SetBinContent(i+1, N)
        N_sim += h_sim.GetBinContent(i+1)
        N_rec += h_rec.GetBinContent(i+1)
        eff = h_eff.GetBinContent(i+1)
        N_total += N
    
        A = scipy.integrate.dblquad(ProjectedArea,
                                    bin.GetXMin(), bin.GetXMax(),
                                    lambda x: bin.GetYMin(),
                                    lambda x: bin.GetYMax())[0]/sa

        try:
            flux_and_error = GetFluxAndError(N,  h_rec.GetBinContent(i+1), h_sim.GetBinContent(i+1), A*sa*dT)
            h_flux.SetBinContent(i+1, flux_and_error[0])
            
            #n_g = g_flux.GetN()
            #g_flux.Set(n_g+1)
            #g_flux.SetPoint(n_g, i, flux_and_error[0])
            #g_flux.SetPointError(n_g, 0, flux_and_error[1])

            h_I.Fill(flux_and_error[0]/flux(theta, 1585), flux_and_error[0]/flux_and_error[1])

            #g_flux_v_sec.Set(n_g+1)
            #g_flux_v_sec.SetPoint(n_g, 1/cos(theta), flux_and_error[0])
            #g_flux_v_sec.SetPointError(n_g, 0, flux_and_error[1])

            #ratio = GetFluxRatio(bin.GetYMin(), bin.GetYMax(),
                                 #bin.GetXMin(), bin.GetXMax())
            #n_g = g_flux_scaled.GetN()
            #g_flux_scaled.Set(n_g+1)
            #g_flux_scaled.SetPoint(n_g, i, flux_and_error[0]/ratio)
            #g_flux_scaled.SetPointError(n_g, 0, flux_and_error[1]/ratio)
            
            flux_and_error = GetFluxAndError(N,  h_rec.GetBinContent(i+1), h_sim.GetBinContent(i+1), A*sa*dT*cos(theta))
            h_vflux.SetBinContent(i+1, flux_and_error[0])
            #n_g = g_vflux.GetN()
            #g_vflux.Set(n_g+1)
            #g_vflux.SetPoint(n_g, i, flux_and_error[0])
            #g_vflux.SetPointError(n_g, 0, flux_and_error[1])

            #ratio = GetVFluxRatio(bin.GetYMin(), bin.GetYMax(),
                                 #bin.GetXMin(), bin.GetXMax())
            #n_g = g_vflux_scaled.GetN()
            #g_vflux_scaled.Set(n_g+1)
            #g_vflux_scaled.SetPoint(n_g, i, flux_and_error[0]/ratio)
            #g_vflux_scaled.SetPointError(n_g, 0, flux_and_error[1]/ratio)

            eff_err = sqrt(N_rec)/N_sim
            
            sum_flux += eff * A * sa
            sum_sq_flux_err += (eff_err * A * sa)**2
            sum_vflux += eff * A * sa * pow(cos(theta),1.00)
            sum_sq_vflux_err += (eff_err * A * sa)**2
        except ZeroDivisionError:
            print("N_rec = %f, N_sim = %f, A = %f"%(h_rec.GetBinContent(i+1), h_sim.GetBinContent(i+1), A))
            print("eff = %f, cos(theta) = %f, sa = %f"%
              (eff, cos(theta), sa))

#A = scipy.integrate.dblquad(ProjectedArea, -pi, pi,
#                            lambda x: first_bin_theta,
#                            lambda x: max_theta)[0]

#tree_mc = ROOT.TChain("tree")
#tree_mc.Add("simulations/muonsim_*.root")
#N_sim_alltheta = 2e6
#N_rec_alltheta = tree_mc.GetEntries("fEventHeader.fTaggedAsMuon")

#tree = ROOT.TChain("tree")
#tree.Add("golden_runs/muon_*.root")
#N_alltheta = tree.GetEntries("fEventHeader.fTaggedAsMuon")
#A_alltheta = scipy.integrate.dblquad(ProjectedArea, -pi, pi,
#                                     lambda x: 0,
#                                     lambda x: pi/2)[0]

flux_in_region = N_total/(sum_flux*dT)
flux_in_region_err = flux_in_region*sqrt(1/N_total
                                         + sum_sq_flux_err/(sum_flux)**2)

print("Total muons in region: %i +/- %0.1f, with acceptance: %0.3e +/- %0.3e."%(N_total, sqrt(N_total), sum_flux, sqrt(sum_sq_flux_err)))
print("Total time: %0.4e s"%(dT))

#total_flux = flux_in_region / 4.416e-01
#total_vflux = flux_in_region / 6.001e-01

total_flux = flux_in_region / 4.134e-01
total_flux_err = flux_in_region_err / 4.134e-01
total_flux_errp = flux_in_region / 4.151e-01 - total_flux
total_flux_errm = total_flux - flux_in_region / 4.115e-01

total_vflux = flux_in_region / 5.063e-01
total_vflux_err = flux_in_region_err / 5.063e-01
total_vflux_errp = flux_in_region / 5.102e-01 - total_vflux
total_vflux_errm = total_vflux - flux_in_region / 5.024e-01

print("Flux in region: (%0.2e +/- %0.2e Hz cm^-2 sr^-1"%(flux_in_region,
                                                         flux_in_region_err))
print("  this scales to a total flux of (%0.2e +/- %0.2e stat + %0.2e - %0.2e sys) Hz cm^-2 sr^-1"%(total_flux, total_flux_err, total_flux_errp, total_flux_errm))
print("  and a vertical flux of (%0.2e +/- %0.2e stat + %0.2e - %0.2e sys) Hz cm^-2 sr^-1"%(total_vflux, total_vflux_err, total_vflux_errp, total_vflux_errm))

#print("%0.2f muons in %0.2f s. Efficiency was %0.3f. Acceptance was %0.3f"%
#      (N_total, dT, N_rec/N_sim, A))
#print("Flux was %0.2e s^-1 cm^-2 sr^-1."%(N_total/(dT*A*N_rec/N_sim)))
#print("Flux was %0.2e s^-1 cm^-2 sr^-1."%(sum_flux/(sum_sa*dT)))
#print("Flux was %0.2e s^-1 cm^-2 sr^-1."%(sum_vflux/(sum_sa*dT)) )

c_sim = ROOT.TCanvas("c_sim")
h_sim.Draw("COLZ")

h_sim_1d = CollapseTheta(h_sim)
c_sim_1d = ROOT.TCanvas("c_sim_1d")
h_sim_1d.Draw()

c_rec = ROOT.TCanvas("c_rec")
h_rec.Draw("COLZ")

c_obs = ROOT.TCanvas("c_obs")
h_obs.Draw("COLZ")

c_eff = ROOT.TCanvas("c_eff")
h_eff.Draw("COLZ")

c_num = ROOT.TCanvas("c_num")
h_num.Draw("COLZ")

c_check = ROOT.TCanvas("c_check")
h_check.Draw("COLZ")

c_flux = ROOT.TCanvas("c_flux")
h_flux.Draw("COLZ")

#h_flux_1d = CollapseTheta(h_flux)
#c_flux_1d = ROOT.TCanvas("c_flux_1d")
#h_flux_1d.Draw()

#ROOT.gStyle.SetOptFit(111)
#fit = ROOT.TF1("fit","[0]*pow(cos(x),1.53)*sin(x)*exp(-0.0008*1585*(1/cos(x)-1))", 0.4, 1.3)
#fit.SetParameter(0, 3e-7)
#fit.SetParameter(1, 1.53)
#fit.SetParName(0, "#Phi")
#fit.SetParName(1, "n")
#h_flux_1d.Fit(fit,"R")

h_flux_bin = Collapse(h_flux)
c_flux_bin = ROOT.TCanvas("c_flux_bin")
h_flux_bin.Draw()
print("Mean: %0.2e +/- %0.2e."%(h_flux_bin.GetMean(),
                                h_flux_bin.GetMeanError()))
f_flux = ROOT.TF1("f_flux", "gaus", 1e-7, 6e-7)
h_flux_bin.Fit(f_flux, "RL")


#c_g_flux = ROOT.TCanvas("c_g_flux")
#g_flux.SetTitle("Muon Flux")
#g_flux.GetYaxis().SetTitle("Flux (s^{-1} cm^{-2} sr^{-1})")
##g_flux.Draw("AP")
#f_g_flux = ROOT.TF1("f_g_flux", "pol0")
#g_flux.Fit(f_g_flux, "")

#c_g_flux_v_sec = ROOT.TCanvas("c_g_flux_v_sec")
#g_flux_v_sec.SetTitle("Muon Flux")
#g_flux_v_sec.GetXaxis().SetTitle("sec#theta")
#g_flux_v_sec.GetYaxis().SetTitle("Flux (s^{-1} cm^{-2} sr^{-1})")
#g_flux_v_sec.Draw("AP")

#c_g_flux_scaled = ROOT.TCanvas("c_g_flux_scaled")
#g_flux_scaled.SetTitle("Scaled Muon Flux")
#g_flux_scaled.GetYaxis().SetTitle("Flux (s^{-1} cm^{-2} sr^{-1})")
#g_flux_scaled.Draw("AP")
#f_g_flux_scaled = ROOT.TF1("f_g_flux_scaled", "pol0")
#g_flux_scaled.Fit(f_g_flux, "")

c_vflux = ROOT.TCanvas("c_vflux")
h_vflux.Draw("COLZ")

c_I = ROOT.TCanvas("c_I")
h_I.Draw("")

#h_vflux_1d = CollapseTheta(h_vflux)
#c_vflux_1d = ROOT.TCanvas("c_vflux_1d")
#h_vflux_1d.Draw()

h_vflux_bin = Collapse(h_vflux)
c_vflux_bin = ROOT.TCanvas("c_vflux_bin")
h_vflux_bin.Draw()
print("Mean: %0.2e +/- %0.2e."%(h_vflux_bin.GetMean(),
                                h_vflux_bin.GetMeanError()))
f_vflux = ROOT.TF1("f_vflux", "gaus", 1e-7, 6e-7)
h_vflux_bin.Fit(f_vflux, "RL")

#c_g_vflux = ROOT.TCanvas("c_g_vflux")
#g_vflux.SetTitle("Vertical Muon Flux")
#g_vflux.GetYaxis().SetTitle("Flux (s^{-1} cm^{-2} sr^{-1})")
#g_vflux.Draw("AP")
#f_g_vflux = ROOT.TF1("f_g_vflux", "pol0")
#g_vflux.Fit(f_g_vflux, "")

#c_g_vflux_scaled = ROOT.TCanvas("c_g_vflux_scaled")
#g_vflux_scaled.SetTitle("Scaled Vertical Muon Flux")
#g_vflux_scaled.GetYaxis().SetTitle("Flux (s^{-1} cm^{-2} sr^{-1})")
#g_vflux_scaled.Draw("AP")
#f_g_vflux_scaled = ROOT.TF1("f_g_vflux_scaled", "pol0")
#g_vflux_scaled.Fit(f_g_flux, "")

save = raw_input("Save? ")
if 'y' in save:
    suffix = raw_input("Suffix? ").rstrip()
    f = ROOT.TFile("rate_histograms_%s.root"%(suffix),"RECREATE")

    h_sim.Write()
    h_sim_1d.Write()
    h_rec.Write()
    h_obs.Write()
    h_eff.Write()
    h_num.Write()
    h_flux.Write()
    h_flux_1d.Write()
    h_vflux.Write()
    h_vflux_1d.Write()
    h_vflux_bin.Write()

    f.Close()
