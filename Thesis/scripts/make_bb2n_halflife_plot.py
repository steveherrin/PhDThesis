import ROOT
import math

# Label, central value, stat err, sys err, use in global fit
measurements = []
measurements.append(('EXO-200 (2011)', 2.11, 0.04, 0.21, True))
measurements.append(('EXO-200 (2012)', 2.23, 0.017, 0.22, True))
measurements.append(('KamLAND-Zen (2012a)', 2.38, 0.02, 0.14, True))
measurements.append(('KamLAND-Zen (2012b)', 2.30, 0.02, 0.11, True))
measurements.append(('EXO-200 (this work)', 1.983, 0.015, 0.045, False))
measurements.reverse()

texts = []
nll_coeffs = [0, 0, 0]

g = ROOT.TGraphErrors(len(measurements))
for (i, m) in enumerate(measurements):
    mean = m[1]
    sigma = m[2]+m[3]
    g.SetPoint(i, mean, i)
    g.SetPointError(i, sigma, 0)
    texts.append(ROOT.TPaveText(1, i - 0.5, 1.8, i + 0.5))
    texts[-1].SetFillColor(ROOT.kWhite)
    texts[-1].SetBorderSize(0)
    texts[-1].AddText(m[0])
    texts[-1].SetTextSize(0.04)
    texts[-1].SetTextAlign(32)

    # We can combine the measurements analytically
    # using a maximum likelihood method. nll of a
    # Gaussian is quadratic
    if m[4]:
        nll_coeffs[0] += 1.0/(2*sigma**2)
        nll_coeffs[1] += -mean/sigma**2
        nll_coeffs[2] += mean**2/(2*sigma**2)

# Minimum of a quadratic
min_x = -nll_coeffs[1]/(2*nll_coeffs[0])
min_nll = nll_coeffs[0]*min_x**2 + nll_coeffs[1]*min_x + nll_coeffs[2]
# Solve for when the value increases by 0.5
nll_coeffs[2] -= min_nll + 0.5
sigma_x = math.sqrt(nll_coeffs[1]**2 - 4*nll_coeffs[0]*nll_coeffs[2])/(2*nll_coeffs[0])
print("Best fit T_{1/2} = (%0.3f +/- %0.3f)x10^{21} yr"%(min_x, sigma_x))

g_bf = ROOT.TGraph(2)
g_bf.SetPoint(0, min_x, -1)
g_bf.SetPoint(1, min_x, len(measurements))
g_bf.SetLineStyle(2)
g_bf.SetLineColor(ROOT.kBlue)

g.SetMarkerStyle(ROOT.kFullCircle)
g.SetTitle("")
g.SetMinimum(-1)
g.SetMaximum(len(measurements))

c = ROOT.TCanvas("c", "c")
c.SetRightMargin(0.02)
c.SetLeftMargin(0.02)
c.SetTopMargin(0.02)
g.Draw("AP")
g.GetXaxis().SetTitle("2#nu#beta#beta T_{1/2} (#times10^{21} yr)")
g.GetXaxis().SetLimits(0.8, 3)
g.GetYaxis().SetTitle("")
g.GetYaxis().SetNdivisions(0, False)
g_bf.Draw("L")

for text in texts:
    text.Draw()

dummy = raw_input("...")
