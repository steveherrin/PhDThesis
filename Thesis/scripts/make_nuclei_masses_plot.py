import ROOT
from math import sqrt, pow

def arrow_label(arrow, text, align):
    x = 0.5*(arrow.GetX1() + arrow.GetX2())
    y = 0.5*(arrow.GetY1() + arrow.GetY2())
    dx = 1.0/20*len(text)
    dy = 1
    dx *= 2-(align/10)
    dy *= (align%10)-2
    t = ROOT.TPaveText(x, y, x+dx, y+dy, "")
    t.SetFillColor(ROOT.kWhite)
    t.SetFillStyle(4000)
    t.SetBorderSize(0)
    t.AddText(text)
    t.SetTextSize(0.03)
    t.SetTextAlign(align)
    return t
    
    

Z0 = 53 # start with iodine 136
M = [135.91465, #136I
     135.907219, #136Xe
     135.9073116, #136Cs
     135.9045759, #136Ba
     135.90764, #136La
     135.907172, #136Ce
     135.912692] #136Pr
label = ["^{136}I ",
         "^{136}Xe",
         "^{136}Cs",
         "^{136}Ba",
         "^{136}La",
         "^{136}Ce",
         "^{136}Pr"]
dE = []

u_to_MeV = 931.454 # MeV/u

E0 = min(M)*u_to_MeV

graphs = []
texts = []
mg = ROOT.TMultiGraph("mg","")

dZ = 0.25 # half width of level lines

for i in xrange(len(M)):
    Z = Z0 + i
    E = M[i]*u_to_MeV - E0
    dE.append(E)

    g = ROOT.TGraph(2)
    g.SetPoint(0, Z-dZ, E)
    g.SetPoint(1, Z+dZ, E)

    g.SetLineWidth(6)
    if Z%2 == 0:
        g.SetLineColor(ROOT.kBlue)
    else:
        g.SetLineColor(ROOT.kGreen)
        
    mg.Add(g, "l")
    graphs.append(g) # So Python doesn't delete them

    t = ROOT.TPaveText(Z-dZ, E-1, Z+dZ, E, "")
    t.AddText(label[i])
    t.SetBorderSize(0)
    t.SetFillColor(ROOT.kWhite)
    t.SetFillStyle(4000)
    t.SetTextSize(0.04)
    texts.append(t)

c = ROOT.TCanvas("c")
c.SetTopMargin(0.02)
c.SetRightMargin(0.02)
mg.Draw("AP")
mg.GetXaxis().SetTitle("Z")
mg.GetYaxis().SetTitle("#DeltaE (MeV)")
mg.GetYaxis().SetRangeUser(-1, 10)

for t in texts:
    t.Draw()

arrows = []
ar_labels = []

ar_sz = 0.01
ar_style = "|>"
dy = 0.1
arrows.append(ROOT.TArrow(53+dZ, dE[0], 54, dE[1]+dy, ar_sz, ar_style))
ar_labels.append(arrow_label(arrows[-1], "#beta^{-}", 13))
arrows.append(ROOT.TArrow(55+dZ, dE[2], 56-1*dZ/3, dE[3]+dy, ar_sz, ar_style))
ar_labels.append(arrow_label(arrows[-1], "#beta^{-}", 13))
arrows.append(ROOT.TArrow(57-dZ, dE[4], 56+1*dZ/3, dE[3]+dy, ar_sz, ar_style))
ar_labels.append(arrow_label(arrows[-1], "#beta^{+}/EC", 33))
arrows.append(ROOT.TArrow(59-dZ, dE[6], 58, dE[5]+dy, ar_sz, ar_style))
ar_labels.append(arrow_label(arrows[-1], "#beta^{+}/EC", 33))

arrows.append(ROOT.TArrow(54+dZ, dE[1], 56-2*dZ/3, dE[3]+dy, ar_sz, ar_style))
ar_labels.append(arrow_label(arrows[-1], "#beta^{-}#beta^{-}", 31))
arrows[-1].SetLineColor(ROOT.kRed)
arrows[-1].SetFillColor(ROOT.kRed)
ar_labels[-1].SetTextColor(ROOT.kRed)
arrows.append(ROOT.TArrow(58-dZ, dE[5], 56+2*dZ/3, dE[3]+dy, ar_sz, ar_style))
ar_labels.append(arrow_label(arrows[-1], "#beta^{+}#beta^{+}/DEC (?)", 11))
arrows[-1].SetLineColor(ROOT.kMagenta)
arrows[-1].SetFillColor(ROOT.kMagenta)
ar_labels[-1].SetTextColor(ROOT.kMagenta)

for ar_label in ar_labels:
    ar_label.Draw()

for arrow in arrows:
    arrow.Draw()

dummy = raw_input("Press Enter...")
