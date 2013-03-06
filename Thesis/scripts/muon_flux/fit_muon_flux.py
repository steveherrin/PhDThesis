#!/bin/env python

import ROOT
ROOT.gSystem.Load("libEXOUtilities")

run_list = open("fileListMasked.dat","r")

tree = ROOT.TChain("tree")

for run_file in run_list.read().splitlines():
    dummy = tree.Add(run_file)

print("Got %i events."%(tree.GetEntries()))

h = ROOT.TH1D("h","Muon Distribution;Zenith Angle;Number/(#pi/100)",
              50, 0, ROOT.TMath.Pi()/2)
#h.Sumw2()
tree.Draw("fEventHeader.fMuonTheta>>h","fEventHeader.fTaggedAsMuon")

dummy = raw_input("Press Enter...")

