import ROOT
from math import pi, sqrt, pow, exp
import scipy.integrate
import numpy
from array import array

alpha = 7.2973e-3
m_e = 0.51099892
Z_Xe = 54
Q = 2.4578

def F(Z, KE):

    E = KE + m_e
    W = E/m_e

    Z0 = Z + 2

    if W <= 1:
        W = 1 + 1e-4

    if W > 2.2:
        a = -8.46e-2 + 2.48e-2*Z0 + 2.37e-4*Z0**2
        b = 1.15e-2 + 3.58e-4*Z0 - 6.17e-5*Z0**2
    else:
        a = -0.811 + 4.46e-2*Z0 + 1.08e-4*Z0**2
        b = 0.673 - 1.82e-2*Z0 + 6.38e-5*Z0**2

    x = sqrt(W-1)
    p = sqrt(W**2 - 1)

    if (p <= 0):
        result = 1
    else:
        result = W/p*exp(a + b*x)

    return result

def D(D, K, i):
    Z = Z_Xe

    T0 = Q/m_e

    E1 = 0.5*(K+D) + 1
    E2 = 0.5*(K+D) + 1

    p1 = sqrt(E1**2 - 1)
    p2 = sqrt(E2**2 - 1)

    T1 = E1 - 1
    T2 = E2 - 1

    return p1*E1*F(Z, T1*m_e)*p2*E2*F(Z, T1*m_e)*pow(T0 - K, i)

def SumSpectrum(K, i):
    if K < 0:
        return 0
    elif K > Q:
        return 0
    a = -K/m_e
    b = K/m_e

    x = scipy.integrate.quad(D, a, b, (K/m_e, i))[0]
    if x < 0:
        return 0
    else:
        return x

def gauss_conv(x, y, res):
    N = len(x)
    
    mu = numpy.mean(x)
    s = res*mu

    gauss = [1.0/(s*sqrt(2*pi))*exp(-0.5*((a-mu)/s)**2) for a in x]

    convolution = numpy.convolve(y, gauss,'same')

    return convolution

def normalize(y, eps, f):
    return [a*f for a in y]

N = 1000
min_E = 0.0
max_E = 1.2
E_scaled = array('d', numpy.linspace(min_E, max_E, N, False))
Es = array('d', (E*Q for E in E_scaled))

eps = (max_E - min_E)/N

bb0n = [0.5/eps if abs(E-Q)<eps else 0 for E in Es]

bb2n = [SumSpectrum(E, 5) for E in Es]

bb0n_smeared = gauss_conv(Es, bb0n, 0.02)
bb2n_smeared = gauss_conv(Es, bb2n, 0.02)

bb0n_int = scipy.integrate.simps(bb0n_smeared, None, eps)
bb0n_norm = array('d', normalize(bb0n_smeared, eps, 1e-2/bb0n_int))
bb2n_int = scipy.integrate.simps(bb2n_smeared, None, eps)
bb2n_norm = array('d', normalize(bb2n_smeared, eps, 1/bb2n_int))

g_bb0n = ROOT.TGraph(N, E_scaled, bb0n_norm)
g_bb0n.SetTitle("")
g_bb0n.SetLineStyle(ROOT.kDashed)
g_bb2n = ROOT.TGraph(N, E_scaled, bb2n_norm)
g_bb2n.SetTitle("")

bb0nX = []
bb0nX.append([0.5/eps if abs(E-Q)<eps else 0 for E in Es])
for i in [1, 2, 3, 5, 7]:
    bb0nX.append([SumSpectrum(E, i) for E in Es])

bb0nX_graphs = []
for bb0nXn in bb0nX:
    bb0nX_int = scipy.integrate.simps(bb0nXn, None, eps)
    bb0nX_norm = array('d', normalize(bb0nXn, eps, 1/bb0nX_int))
    g_bb0nX = ROOT.TGraph(N, E_scaled, bb0nX_norm)
    bb0nX_graphs.append(g_bb0nX)

min_E = 0.9
max_E = 1.1
E_scaled_z = array('d', numpy.linspace(min_E, max_E, N, False))
Es_z = array('d', (E*Q for E in E_scaled_z))

eps_z = (max_E - min_E)/N

bb0n_z = [0.5/eps_z if abs(E-Q)<eps_z else 0 for E in Es_z]

bb2n_z = [SumSpectrum(E, 5) for E in Es_z]

bb0n_smeared_z = gauss_conv(Es_z, bb0n_z, 0.02)
bb2n_smeared_z = gauss_conv(Es_z, bb2n_z, 0.02)

bb0n_norm_z = array('d', normalize(bb0n_smeared_z, eps, 1e-6/bb0n_int))
bb2n_norm_z = array('d', normalize(bb2n_smeared_z, eps, 1.0/bb2n_int))

g_bb0n_z = ROOT.TGraph(N, E_scaled_z, bb0n_norm_z)
g_bb0n_z.SetTitle("")
g_bb0n_z.SetLineStyle(ROOT.kDashed)
g_bb2n_z = ROOT.TGraph(N, E_scaled_z, bb2n_norm_z)
g_bb2n_z.SetTitle("")

#print("bb0n %f"%(sum((y*eps for y in bb0n_norm))))
#print("bb2n %f"%(sum((y*eps for y in bb2n_norm))))

c_both = ROOT.TCanvas("c_both","c_both")
p = ROOT.TPad("p", "p", 0, 0, 1, 1)
p.SetRightMargin(0.02)
p.SetTopMargin(0.02)
p.Draw()
p.cd()
g_bb2n.Draw("AL")
g_bb0n.Draw("L")
g_bb2n.GetYaxis().SetTitle("dN/dE")
g_bb2n.GetXaxis().SetTitle("Sum e^{-} Energy (E/Q)")

c_both.cd()
p_inset = ROOT.TPad("p_inset","p_inset",0.5, 0.5, 0.995, 0.995)
p_inset.SetRightMargin(0.05)
p_inset.SetTopMargin(0.05)
p_inset.Draw()
p_inset.cd()
g_bb2n_z.Draw("AL")
g_bb0n_z.Draw("L")
g_bb2n_z.GetYaxis().SetTitle("dN/dE")
g_bb2n_z.GetXaxis().SetTitle("Sum e^{-} Energy (E/Q)")
g_bb2n_z.GetYaxis().SetNoExponent(False)
# Zoom in so we can't see edge effects of the convolution
g_bb2n_z.GetXaxis().SetRangeUser(1-0.25*(1-min_E), 1+0.25*(max_E-1))
g_bb2n_z.GetYaxis().SetRangeUser(0, 0.0004)

c_z = ROOT.TCanvas("c_z","c_z")
c_z.SetRightMargin(0.05)
c_z.SetTopMargin(0.05)
g_bb2n_z.Draw("AL")
g_bb0n_z.Draw("L")

c = ROOT.TCanvas("c","c")
c.SetRightMargin(0.05)
c.SetTopMargin(0.05)
g_bb2n.Draw("AL")
g_bb0n.Draw("L")

c_majoron = ROOT.TCanvas("c_majoron")
c_majoron.SetRightMargin(0.05)
c_majoron.SetTopMargin(0.05)
colors = [ROOT.kBlack, ROOT.kRed, ROOT.kGreen, ROOT.kBlue,
          ROOT.kMagenta, ROOT.kCyan]

draw_opt = "AL"
for i in xrange(len(bb0nX_graphs)):
    bb0nX_graphs[-(i+1)].SetLineColor(colors[-(i+1)])
    bb0nX_graphs[-(i+1)].Draw(draw_opt)
    draw_opt = "L"
# Draw bb0n last so it doesn't scale others to 0
bb0nX_graphs[-1].SetTitle("")
bb0nX_graphs[-1].GetXaxis().SetRangeUser(0, 1.1)
bb0nX_graphs[-1].GetXaxis().SetTitle("Sum e^{-} Energy (E/Q)")
bb0nX_graphs[-1].GetYaxis().SetTitle("dN/dE")

l_majoron = ROOT.TLegend(0.45, 0.77, 0.85, 0.94)
l_majoron.SetFillColor(ROOT.kWhite)
l_majoron.SetNColumns(2)
l_majoron.AddEntry(bb0nX_graphs[0], "0#nu#beta#beta", "l")
l_majoron.AddEntry(bb0nX_graphs[1], "0#nu#beta#beta#chi^{0} (n=1)", "l")
l_majoron.AddEntry(bb0nX_graphs[4], "2#nu#beta#beta (n=5)", "l")
l_majoron.AddEntry(bb0nX_graphs[2], "0#nu#beta#beta#chi^{0} (n=2)", "l")
l_majoron.AddEntry(None, "", "")
l_majoron.AddEntry(bb0nX_graphs[3], "0#nu#beta#beta#chi^{0}(#chi^{0}) (n=3)", "l")
l_majoron.AddEntry(None, "", "")
l_majoron.AddEntry(bb0nX_graphs[5], "0#nu#beta#beta#chi^{0}#chi^{0} (n=7)", "l")
l_majoron.Draw()

dummy = raw_input("Press Enter...")
