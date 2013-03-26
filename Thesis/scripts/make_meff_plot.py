import ROOT
import numpy
from math import sqrt, fabs, log10

class Parameter:
    def __init__(self, bf, m1s, p1s, m2s, p2s, m3s, p3s):
        self.bf = bf
        self.m = [bf, m1s, m2s, m3s]
        self.p = [bf, p1s, p2s, p3s]
        self.CalcDeltas()
    def CalcDeltas(self):
        self.d_m = [self.bf - m for m in self.m]
        self.d_p = [p - self.bf for p in self.p]

# From arXiv:1205.4018
# Parameters are best fit, then 1, 2, 3 sigma RANGES
# If you want to use PLUS/MINUS, you need to modify
# the class above (or better, introduce a new one)
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

# limits
exo_limit = 0.140

# 0.23 eV >> mass splittings, so could expand to zeroth order and call it a day
cosmo_sum = 0.23
#cosmo_limit = cosmo_sum/3
# but let's use the 1st order expansion, which Mathematica tells me is
cosmo_limit_ih = cosmo_sum/3 - m_atm_ih.bf/cosmo_sum - 0.5*m_sol.bf/cosmo_sum
cosmo_limit_nh = cosmo_sum/3 - 0.5/cosmo_sum*(m_atm_nh.bf + m_sol.bf)
cosmo_limit = max(cosmo_limit_ih, cosmo_limit_nh)

#print("IH cosmic limit 0th: %0.3e 1st: %0.3e"%(cosmo_limit, cosmo_limit_ih))
#print("NH cosmic limit 0th: %0.3e 1st: %0.3e"%(cosmo_limit, cosmo_limit_nh))      


# plot range
log_m_lo = log10(0.9e-4)
log_m_hi = log10(1.1)
n = 2000

graph_lo = 1e-4
graph_hi = 1e0

def m_sum2(m1, m2):
    return sqrt(m1**2 + m2)

def m_sum3(m1, m2, m3):
    return sqrt(m1**2 + m2 + m3)

def sigma_m_sol_nh(m, s_12, s_23, s_13, m_sol, m_atm, pm_1, pm_2, d_m_sol):
    d = pm_1*(1-s_12)*(1-s_13)/(2*m_sum2(m, m_sol))
    return d * d_m_sol

def sigma_m_sol_ih(m, s_12, s_23, s_13, m_sol, m_atm, pm_1, pm_2, d_m_sol):
    d = (1-s_12**2)*(1-s_13)/(2*m_sum3(m, m_sol, m_atm))
    return d * d_m_sol

def sigma_s_12_nh(m, s_12, s_23, s_13, m_sol, m_atm, pm_1, pm_2, d_s_12):
    d = (-m + pm_1*m_sum2(m, m_atm))*(1-s_13)
    return d * d_s_12

def sigma_s_12_ih(m, s_12, s_23, s_13, m_sol, m_atm, pm_1, pm_2, d_s_12):
    d = (-m_sum2(m, m_atm) + pm_1*m_sum3(m, m_atm, m_sol))*(1-s_13)
    return d * d_s_12

def sigma_m_atm_nh(m, s_12, s_23, s_13, m_sol, m_atm, pm_1, pm_2, d_m_atm):
    d = pm_2*s_13/(2*m_sum2(m, m_atm))
    return d * d_m_atm

def sigma_m_atm_ih(m, s_12, s_23, s_13, m_sol, m_atm, pm_1, pm_2, d_m_atm):
    d = (1-s_12**2)*(1-s_13)/(2*m_sum2(m, m_atm))
    return d * d_m_atm

def sigma_s_13_nh(m, s_12, s_23, s_13, m_sol, m_atm, pm_1, pm_2, d_s_13):
    d = -m*(1-s_12)-pm_1*m_sum2(m, m_sol)+pm_2*m_sum2(m, m_atm)
    return d * d_s_13

def sigma_s_13_ih(m, s_12, s_23, s_13, m_sol, m_atm, pm_1, pm_2, d_s_13):
    d = -m_sum2(m, m_atm)*(1-s_12)-pm_1*m_sum3(m, m_atm, m_sol)+pm_2*m
    return d * d_s_13

def get_ih(m, pm_1, pm_2, err_dir):

    result = [0, 0, 0, 0]
    
    result[0] = (m_sum2(m, m_atm_ih.bf)*(1-s_12.bf)*(1-s_13_ih.bf)
                 + pm_1 * m_sum3(m, m_atm_ih.bf, m_sol.bf)*s_12.bf*(1-s_13_ih.bf)
                 + pm_2 * m * s_13_ih.bf)
    if result[0] < 0:
        err_dir *= -1

    for n_sigma in xrange(1, 3+1):

        sum_sq_sigma = 0

        for (f, x) in [(sigma_m_sol_ih, m_sol),
                       (sigma_m_atm_ih, m_atm_ih),
                       (sigma_s_12_ih, s_12),
                       (sigma_s_13_ih, s_13_ih)]:

            added_sigma = 0
            for dx in [x.d_m[n_sigma], x.d_p[n_sigma]]:
    
                sigma = f(m, s_12.bf, s_23_ih.bf,
                          s_13_ih.bf, m_sol.bf,
                          m_atm_ih.bf, pm_1, pm_2, dx)
                
                if err_dir > 0 and sigma > added_sigma:
                    added_sigma = sigma
                elif err_dir < 0 and sigma < added_sigma:
                    added_sigma = sigma

            sum_sq_sigma += added_sigma**2

        if err_dir > 0:
            result[n_sigma] = result[0] + sqrt(sum_sq_sigma)
        else:
            result[n_sigma] = result[0] - sqrt(sum_sq_sigma)
      
    return result


def get_nh(m, pm_1, pm_2, err_dir):

    result = [0, 0, 0, 0]
    
    result[0] = (m*(1-s_12.bf)*(1-s_13_nh.bf)
                 + pm_1 * m_sum2(m, m_sol.bf)*s_12.bf*(1-s_13_nh.bf)
                 + pm_2 * m_sum2(m, m_atm_nh.bf) * s_13_nh.bf)


    for n_sigma in xrange(1, 3+1):

        sum_sq_sigma = 0

        for (f, x) in [(sigma_m_sol_nh, m_sol),
                       (sigma_m_atm_nh, m_atm_nh),
                       (sigma_s_12_nh, s_12),
                       (sigma_s_13_nh, s_13_nh)]:

            added_sigma = 0
            for dx in [x.d_m[n_sigma], x.d_p[n_sigma]]:
    
                sigma = f(m, s_12.bf, s_23_nh.bf,
                          s_13_nh.bf, m_sol.bf,
                          m_atm_nh.bf, pm_1, pm_2, dx)

                if err_dir > 0 and sigma > added_sigma:
                    added_sigma = sigma
                elif err_dir < 0 and sigma < added_sigma:
                    added_sigma = sigma
            sum_sq_sigma += added_sigma**2
        if err_dir > 0:
            result[n_sigma] = result[0] + sqrt(sum_sq_sigma)
        else:
            result[n_sigma] = result[0] - sqrt(sum_sq_sigma)
      
    return result

m_l = []



nh_lo_l = []
nh_hi_l = []

ih_lo_l = []
ih_hi_l = []

for i, m in enumerate(numpy.logspace(log_m_lo, log_m_hi, n, False)):
    nh_lo_mm = get_nh(m, -1, -1, -1)
    # yes, positive error on this one because this curve only
    # shows up when the expression is negative
    nh_lo_mp = get_nh(m, -1, 1, 1)

    nh_lo = []
    for (y1, y2) in zip(nh_lo_mm, nh_lo_mp):
        if y1 > 0:
            nh_lo.append(y1)
        elif y1 < 0 and y2 < 0:
            # colvoluted, and there's probably a better way
            nh_lo.append(fabs(nh_lo_mp[0])-fabs(y2-nh_lo_mp[0]))
        else:
            nh_lo.append(0)
    
    nh_hi_pp = get_nh(m, 1, 1, 1)
    nh_hi_pm = get_nh(m, 1, -1, 1)
    nh_hi = [max(y1, y2) for (y1, y2) in zip(nh_hi_pp, nh_hi_pm)]
    
    ih_lo_mm = get_ih(m, -1, -1, -1)
    ih_lo_mp = get_ih(m, -1, 1, -1)
    ih_lo = [min(y1, y2) for (y1, y2) in zip(ih_lo_mm, ih_lo_mp)]
    ih_hi_pp = get_ih(m, 1, 1, 1)
    ih_hi_pm = get_ih(m, 1, -1, 1)
    ih_hi = [max(y1, y2) for (y1, y2) in zip(ih_hi_pp, ih_hi_pm)]

    m_l.append(m)
    nh_lo_l.append(nh_lo)
    nh_hi_l.append(nh_hi)
    ih_lo_l.append(ih_lo)
    ih_hi_l.append(ih_hi)

print("For a limit on m_bb of %0.3e eV:"%(exo_limit))
found_limit = [False for x in xrange(4)]
for (i, m) in enumerate(m_l):
    for j in xrange(1, len(found_limit)):
        y = min(nh_lo_l[i][j], ih_lo_l[i][j])
        if (not found_limit[j]) and y > exo_limit:
            print("    %i sigma limit on m_min: %0.3e eV"%(j, m))
            found_limit[j] = True
            
    
    
g_nh_1s = ROOT.TGraph(2*n)
g_nh_2s = ROOT.TGraph(2*n)
g_nh_3s = ROOT.TGraph(2*n)
g_ih_1s = ROOT.TGraph(2*n)
g_ih_2s = ROOT.TGraph(2*n)
g_ih_3s = ROOT.TGraph(2*n)

for i in range(n):
    g_nh_1s.SetPoint(i, m_l[i], nh_hi_l[i][1])
    g_nh_1s.SetPoint(n+i, m_l[n-i-1], nh_lo_l[n-i-1][1])
    #g_nh_1s.SetPoint(n+i, m_l[n-i-1], 0)
    g_nh_2s.SetPoint(i, m_l[i], nh_hi_l[i][2])
    g_nh_2s.SetPoint(n+i, m_l[n-i-1], nh_lo_l[n-i-1][2])
    #g_nh_2s.SetPoint(n+i, m_l[n-i-1], 0)
    g_nh_3s.SetPoint(i, m_l[i], nh_hi_l[i][3])
    g_nh_3s.SetPoint(n+i, m_l[n-i-1], nh_lo_l[n-i-1][3])
    #g_nh_3s.SetPoint(n+i, m_l[n-i-1], 0)
    g_ih_1s.SetPoint(i, m_l[i], ih_hi_l[i][1])
    g_ih_1s.SetPoint(n+i, m_l[n-i-1], ih_lo_l[n-i-1][1])
    g_ih_2s.SetPoint(i, m_l[i], ih_hi_l[i][2])
    g_ih_2s.SetPoint(n+i, m_l[n-i-1], ih_lo_l[n-i-1][2])
    g_ih_3s.SetPoint(i, m_l[i], ih_hi_l[i][3])
    g_ih_3s.SetPoint(n+i, m_l[n-i-1], ih_lo_l[n-i-1][3])

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
c.SetRightMargin(0.02)
c.SetTopMargin(0.02)

g_nh_3s.Draw("AF")
g_nh_2s.Draw("F")
g_nh_1s.Draw("F")
g_ih_3s.Draw("F")
g_ih_2s.Draw("F")
g_ih_1s.Draw("F")

l = ROOT.TLegend(.85, .15, .95, .30)
l.SetFillColor(ROOT.kWhite)
l.AddEntry(g_nh_1s, "1 #sigma", "f")
l.AddEntry(g_nh_2s, "2 #sigma", "f")
l.AddEntry(g_nh_3s, "3 #sigma", "f")
l.Draw()

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
g_exo.SetPoint(0, graph_lo, exo_limit)
g_exo.SetPoint(1, graph_hi, exo_limit)
g_exo.SetLineColor(ROOT.kRed)
g_exo.SetLineStyle(2)
g_exo.Draw()
t_exo = ROOT.TText(6e-4, .16, "EXO-200 Limit")
t_exo.SetTextColor(ROOT.kRed)
t_exo.SetTextSize(0.05)
t_exo.SetTextFont(42)
t_exo.Draw()

#print("Cosmo limit was %f."%(cosmo_limit))

g_cosmo = ROOT.TGraph(2)
g_cosmo.SetPoint(0, cosmo_limit, graph_lo)
g_cosmo.SetPoint(1, cosmo_limit, graph_hi)
g_cosmo.SetLineColor(ROOT.kRed)
g_cosmo.SetLineStyle(2)
g_cosmo.Draw()
t_cosmo = ROOT.TText(cosmo_limit+0.02, 2e-2, "Cosmological Limit")
t_cosmo.SetTextColor(ROOT.kRed)
t_cosmo.SetTextSize(0.05)
t_cosmo.SetTextFont(42)
t_cosmo.SetTextAngle(-90)
t_cosmo.Draw()

raw_input("Press enter...")
