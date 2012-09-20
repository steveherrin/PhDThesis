import ROOT
import numpy
from math import sqrt, fabs, log10

class Parameter:
    def __init__(self, bf, m1s, p1s, m2s, p2s, m3s, p3s):
        self.bf = bf
        self.m1s = m1s
        self.p1s = p1s
        self.m2s = m2s
        self.p2s = p2s
        self.m3s = m3s
        self.p3s = p3s

class Result:
    def __init__(self):
        self.s1 = 0
        self.s2 = 0
        self.s3 = 0

# From arXiv:1205.4018
m_sol = Parameter(7.62e-5, 7.43e-5, 7.81e-5, 7.27e-5,
                  8.01e-5, 7.12e-5, 8.20e-5)
m_atm_nh = Parameter(2.55e-3, 2.46e-3, 2.61e-3, 2.38e-3,
                     2.68e-3, 2.31e-3, 2.74e-3)
m_atm_ih = Parameter(2.43e-3, 2.37e-3, 2.50e-3, 2.29e-3,
                     2.58e-3, 2.21e-3, 2.64e-3)
s_12 = Parameter(0.320, 0.303, 0.336, 0.290, 0.350, 0.270, 0.370)
s_23_nh = Parameter(0.613, 0.400, 0.635, 0.380, 0.660, 0.360, 0.680)
s_23_ih = Parameter(0.600, 0.569, 0.626, 0.390, 0.650, 0.370, 0.670)
s_13_nh = Parameter(0.0246, 0.0218, 0.0275, 0.019, 0.030, 0.017, 0.033)
s_13_ih = Parameter(0.0250, 0.0223, 0.0276, 0.020, 0.030, 0.017, 0.033)



# plot range
log_m_lo = log10(0.9e-4)
log_m_hi = log10(1.1)
n = 1000

def get_ih_hi(m):
    result = Result()
    result.s1 = (m*s_13_ih.m1s
                 + sqrt(m**2 + m_atm_ih.p1s)*s_12.p1s
                 + sqrt(m**2 + m_atm_ih.p1s - m_sol.m1s)*(1-s_12.p1s)*(1-s_13_ih.m1s))
    result.s2 = (m*s_13_ih.m2s
                 + sqrt(m**2 + m_atm_ih.p2s)*s_12.p2s
                 + sqrt(m**2 + m_atm_ih.p2s - m_sol.m2s)*(1-s_12.p2s)*(1-s_13_ih.m2s))
    result.s3 = (m*s_13_ih.m3s
                 + sqrt(m**2 + m_atm_ih.p3s)*s_12.p3s
                 + sqrt(m**2 + m_atm_ih.p3s - m_sol.m3s)*(1-s_12.p3s)*(1-s_13_ih.m3s))
    return result

def get_ih_lo(m):
    result = Result()
    result.s1 = fabs(m*s_13_ih.p1s
                 + sqrt(m**2 + m_atm_ih.m1s)*s_12.p1s
                 - sqrt(m**2 + m_atm_ih.m1s - m_sol.p1s)*(1-s_12.p1s)*(1-s_13_ih.p1s))
    result.s2 = fabs(m*s_13_ih.p2s
                 + sqrt(m**2 + m_atm_ih.m2s)*s_12.p2s
                 - sqrt(m**2 + m_atm_ih.m2s - m_sol.p2s)*(1-s_12.p2s)*(1-s_13_ih.p2s))
    result.s3 = fabs(m*s_13_ih.p3s
                 + sqrt(m**2 + m_atm_ih.m3s)*s_12.p3s
                 - sqrt(m**2 + m_atm_ih.m3s - m_sol.p3s)*(1-s_12.p3s)*(1-s_13_ih.p3s))
    return result

def get_nh_hi(m):
    result = Result()
    result.s1 = (m*(1-s_12.p1s)*(1-s_13_nh.p1s)
                 + sqrt(m**2 + m_sol.p1s)*s_12.p1s*(1-s_13_nh.p1s)
                 + sqrt(m**2 + m_atm_nh.p1s)*s_13_nh.p1s)
    result.s2 = (m*(1-s_12.p2s)*(1-s_13_nh.p2s)
                 + sqrt(m**2 + m_sol.p2s)*s_12.p2s*(1-s_13_nh.p2s)
                 + sqrt(m**2 + m_atm_nh.p2s)*s_13_nh.p2s)
    result.s3 = (m*(1-s_12.p3s)*(1-s_13_nh.p3s)
                 + sqrt(m**2 + m_sol.p3s)*s_12.p3s*(1-s_13_nh.p3s)
                 + sqrt(m**2 + m_atm_nh.p3s)*s_13_nh.p3s)
    return result

def get_nh_lo(m):
    result = Result()
    result.s1 = (m*(1-s_12.p1s)*(1-s_13_nh.p1s)
                 - sqrt(m**2 + m_sol.p1s)*s_12.p1s*(1-s_13_nh.p1s)
                 - sqrt(m**2 + m_atm_nh.p1s)*s_13_nh.p1s)
    if result.s1 < 0:
        result.s1 = -(m*(1-s_12.m1s)*(1-s_13_nh.p1s)
                      - sqrt(m**2 + m_sol.m1s)*s_12.m1s*(1-s_13_nh.p1s)
                      + sqrt(m**2 + m_atm_nh.p1s)*s_13_nh.p1s)
        if result.s1 < 0:
            result.s1 = 0
    result.s2 = (m*(1-s_12.p2s)*(1-s_13_nh.p2s)
                 - sqrt(m**2 + m_sol.p2s)*s_12.p2s*(1-s_13_nh.p2s)
                 - sqrt(m**2 + m_atm_nh.p2s)*s_13_nh.p2s)
    if result.s2 < 0:
        result.s2 = -(m*(1-s_12.m2s)*(1-s_13_nh.p2s)
                      - sqrt(m**2 + m_sol.m2s)*s_12.m2s*(1-s_13_nh.p2s)
                      + sqrt(m**2 + m_atm_nh.p2s)*s_13_nh.p2s)
        if result.s2 < 0:
            result.s2 = 0
    result.s3 = (m*(1-s_12.p3s)*(1-s_13_nh.p3s)
                 - sqrt(m**2 + m_sol.p3s)*s_12.p3s*(1-s_13_nh.p3s)
                 - sqrt(m**2 + m_atm_nh.p3s)*s_13_nh.p3s)
    if result.s3 < 0:
        result.s3 = -(m*(1-s_12.m3s)*(1-s_13_nh.p3s)
                      - sqrt(m**2 + m_sol.m3s)*s_12.m3s*(1-s_13_nh.p3s)
                      + sqrt(m**2 + m_atm_nh.p3s)*s_13_nh.p3s)
        if result.s3 < 0:
            result.s3 = 0
    return result

m_l = []

nh_lo_l = []
nh_hi_l = []

ih_lo_l = []
ih_hi_l = []

for i, m in enumerate(numpy.logspace(log_m_lo, log_m_hi, n, False)):
    nh_lo = get_nh_lo(m)
    nh_hi = get_nh_hi(m)
    ih_lo = get_ih_lo(m)
    ih_hi = get_ih_hi(m)

    m_l.append(m)
    nh_lo_l.append(nh_lo)
    nh_hi_l.append(nh_hi)
    ih_lo_l.append(ih_lo)
    ih_hi_l.append(ih_hi)
    
g_nh_1s = ROOT.TGraph(2*n)
g_nh_2s = ROOT.TGraph(2*n)
g_nh_3s = ROOT.TGraph(2*n)
g_ih_1s = ROOT.TGraph(2*n)
g_ih_2s = ROOT.TGraph(2*n)
g_ih_3s = ROOT.TGraph(2*n)

for i in range(n):
    g_nh_1s.SetPoint(i, m_l[i], nh_hi_l[i].s1)
    g_nh_1s.SetPoint(n+i, m_l[n-i-1], nh_lo_l[n-i-1].s1)
    #g_nh_1s.SetPoint(n+i, m_l[n-i-1], 0)
    g_nh_2s.SetPoint(i, m_l[i], nh_hi_l[i].s2)
    g_nh_2s.SetPoint(n+i, m_l[n-i-1], nh_lo_l[n-i-1].s2)
    #g_nh_2s.SetPoint(n+i, m_l[n-i-1], 0)
    g_nh_3s.SetPoint(i, m_l[i], nh_hi_l[i].s3)
    g_nh_3s.SetPoint(n+i, m_l[n-i-1], nh_lo_l[n-i-1].s3)
    #g_nh_3s.SetPoint(n+i, m_l[n-i-1], 0)
    g_ih_1s.SetPoint(i, m_l[i], ih_hi_l[i].s1)
    g_ih_1s.SetPoint(n+i, m_l[n-i-1], ih_lo_l[n-i-1].s1)
    g_ih_2s.SetPoint(i, m_l[i], ih_hi_l[i].s2)
    g_ih_2s.SetPoint(n+i, m_l[n-i-1], ih_lo_l[n-i-1].s2)
    g_ih_3s.SetPoint(i, m_l[i], ih_hi_l[i].s3)
    g_ih_3s.SetPoint(n+i, m_l[n-i-1], ih_lo_l[n-i-1].s3)

g_nh_1s.SetFillStyle(1001)
g_nh_1s.SetFillColor(ROOT.kBlue+2)
g_nh_2s.SetFillStyle(1001)
g_nh_2s.SetFillColor(ROOT.kBlue+1)
g_nh_3s.SetFillStyle(1001)
g_nh_3s.SetFillColor(ROOT.kBlue+0)

g_ih_1s.SetFillStyle(1001)
g_ih_1s.SetFillColor(ROOT.kBlue+2)
g_ih_2s.SetFillStyle(1001)
g_ih_2s.SetFillColor(ROOT.kBlue+1)
g_ih_3s.SetFillStyle(1001)
g_ih_3s.SetFillColor(ROOT.kBlue+0)

c = ROOT.TCanvas("c", "Hierarchy", 600, 600)

g_nh_3s.Draw("AF")
g_nh_2s.Draw("F")
g_nh_1s.Draw("F")
g_ih_3s.Draw("F")
g_ih_2s.Draw("F")
g_ih_1s.Draw("F")

l = ROOT.TLegend(.60, .15, .75, .30)
l.SetFillColor(ROOT.kWhite)
l.AddEntry(g_nh_1s, "1 #sigma", "f")
l.AddEntry(g_nh_2s, "2 #sigma", "f")
l.AddEntry(g_nh_3s, "3 #sigma", "f")
l.Draw()

graph_lo = 1e-4
graph_hi = 1e0

c.SetLogy()
c.SetLogx()
g_nh_3s.SetTitle("")
x_axis = g_nh_3s.GetXaxis()
y_axis = g_nh_3s.GetYaxis()
x_axis.SetTitle("m_{min} (eV)")
y_axis.SetTitle("#LT m_{#beta #beta} #GT (eV)")
x_axis.SetTitleSize(0.03)
y_axis.SetTitleSize(0.03)
x_axis.SetTitleOffset(1.5)
y_axis.SetTitleOffset(1.5)
x_axis.SetLimits(graph_lo, graph_hi)
x_axis.SetRangeUser(graph_lo, graph_hi)
y_axis.SetRangeUser(graph_lo, graph_hi)
c.RedrawAxis()

t_nh = ROOT.TText(5e-4, 1.8e-3, "Normal")
t_nh.SetTextColor(ROOT.kWhite)
t_nh.SetTextSize(0.05)
t_nh.SetTextFont(42)
t_nh.Draw()

t_ih = ROOT.TText(5e-4, 2.2e-2, "Inverted")
t_ih.SetTextColor(ROOT.kWhite)
t_ih.SetTextSize(0.05)
t_ih.SetTextFont(42)
t_ih.Draw()

g_exo = ROOT.TGraph(2)
g_exo.SetPoint(0, graph_lo, .140)
g_exo.SetPoint(1, graph_hi, .140)
g_exo.SetLineColor(ROOT.kRed)
g_exo.SetLineStyle(2)
g_exo.Draw()
t_exo = ROOT.TText(8e-4, .16, "EXO-200 Limit")
t_exo.SetTextColor(ROOT.kRed)
t_exo.SetTextSize(0.05)
t_exo.SetTextFont(42)
t_exo.Draw()


g_cosmo = ROOT.TGraph(2)
g_cosmo.SetPoint(0, 0.3, graph_lo)
g_cosmo.SetPoint(1, 0.3, graph_hi)
g_cosmo.SetLineColor(ROOT.kRed)
g_cosmo.SetLineStyle(2)
g_cosmo.Draw()
t_cosmo = ROOT.TText(.35, 3e-2, "Cosmological Limit")
t_cosmo.SetTextColor(ROOT.kRed)
t_cosmo.SetTextSize(0.05)
t_cosmo.SetTextFont(42)
t_cosmo.SetTextAngle(-90)
t_cosmo.Draw()

raw_input("Press enter...")
