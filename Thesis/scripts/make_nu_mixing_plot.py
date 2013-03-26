import ROOT
import numpy
from numpy import dot
from math import sqrt, pi

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
delta_nh = 0.8*pi
delta_ih = -0.03*pi

# Kind of bad notation. Sorry.
s12 = sqrt(s_12.bf)
c12 = sqrt(1-s_12.bf)
s13_nh = sqrt(s_13_nh.bf)
c13_nh = sqrt(1-s_13_nh.bf)
s13_ih = sqrt(s_13_ih.bf)
c13_ih = sqrt(1-s_13_ih.bf)
s23_nh = sqrt(s_23_nh.bf)
c23_nh = sqrt(1-s_23_nh.bf)
s23_ih = sqrt(s_23_ih.bf)
c23_ih = sqrt(1-s_23_ih.bf)

# CP violating phase
cp_nh = numpy.exp(numpy.complex(0,delta_nh))
cp_ih = numpy.exp(numpy.complex(0,delta_ih))

# Factors in the mixing matrix, ignoring
# the Majorana phases because they don't
# matter for oscillation

U12 = numpy.array([[c12, s12, 0],
                   [-s12, c12, 0],
                   [0, 0, 1]])
U13_nh = numpy.array([[c13_nh, 0, s13_nh*cp_nh],
                      [0, 1, 0],
                      [-s13_nh*cp_nh, 0, c13_nh]])
U13_ih = numpy.array([[c13_ih, 0, s13_ih*cp_ih],
                      [0, 1, 0],
                      [-s13_ih*cp_ih, 0, c13_ih]])
U23_nh = numpy.array([[1, 0, 0],
                      [0, c23_nh, s23_nh],
                      [0, -s23_nh, c23_nh]])
U23_ih = numpy.array([[1, 0, 0],
                      [0, c23_ih, s23_ih],
                      [0, -s23_ih, c23_ih]])

# Multiply out the matrices

U_nh = dot(U23_nh, dot(U13_nh, U12))
U_ih = dot(U23_ih, dot(U13_ih, U12))

# o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-
# Draw the colored lines for flavor content

colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen]
line_px = 10

y_min = 0.25
dy21 = 0.1
dy31 = 0.4
ys_nh = [y_min, y_min + dy21, y_min + dy31]
ys_ih = [y_min + dy31, y_min + dy31 + dy21, y_min]

lines = []

for i in xrange(3):

    y_nh = ys_nh[i]
    y_ih = ys_ih[i]

    x_nh = 2.0/12
    x_ih = 7.0/12
    scale = 3.0/12

    for j in xrange(3):
        
        # flavor fraction
        p_nh = scale*abs(U_nh[j, i])**2
        p_ih = scale*abs(U_ih[j, i])**2

        l_nh = ROOT.TLine(x_nh, y_nh,
                          x_nh + p_nh, y_nh)
        
        l_nh.SetLineColor(colors[j])
        l_nh.SetLineWidth(line_px)
        x_nh += p_nh
        lines.append(l_nh)

        l_ih = ROOT.TLine(x_ih, y_ih,
                          x_ih + p_ih, y_ih)
        
        l_ih.SetLineColor(colors[j])
        l_ih.SetLineWidth(line_px)
        x_ih += p_ih
        lines.append(l_ih)

x_size = 640
y_size = 480
c = ROOT.TCanvas("c", "c", x_size, y_size)

for line in lines:
    line.Draw()

# o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-
# Add a legend

legend = ROOT.TLegend(11.0/24, 0.8, 13.0/24, 0.97)
legend.SetFillColor(ROOT.kWhite)
legend.AddEntry(lines[0], "#nu_{e}", "L")
legend.AddEntry(lines[2], "#nu_{#mu}", "L")
legend.AddEntry(lines[4], "#nu_{#tau}", "L")
legend.SetTextAlign(22)
legend.Draw()

# o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-
# Now draw arrows between mass states

x_nh = 2.0/12
x_ih = 7.0/12

ar_size = 0.015

arrows = []
splittings = []
ar_x = 5.0/24

ax_y = 0.05

line_w = 0.5*line_px/x_size

# So the arrows are staggered nicely
ar_ys_nh = [(ax_y-line_w, ys_nh[0]),
            (ys_nh[0], ys_nh[2]),
            (ys_nh[0], ys_nh[1])]

ar_ys_ih = [(ax_y-line_w, ys_ih[2]),
            (ys_ih[2], ys_ih[0]),
            (ys_ih[0], ys_ih[1])]

for i in xrange(3):

    # A little space so they don't overlap
    dx = (i+1.0)/48
        
    ar_nh = ROOT.TArrow(x_nh + dx, ar_ys_nh[i][0]+line_w,
                        x_nh + dx, ar_ys_nh[i][1]-line_w,
                        ar_size, "<|>")
    ar_ih = ROOT.TArrow(x_ih + dx, ar_ys_ih[i][0]+line_w,
                        x_ih + dx, ar_ys_ih[i][1]-line_w,
                        ar_size, "<|>")

    arrows.append(ar_nh)
    arrows.append(ar_ih)

for arrow in arrows:
    arrow.Draw()

# o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-
# Label those arrows

splitting_txts_nh = ["?", "#Deltam_{21}^{2} (sol)", "#Deltam_{31}^{2} (atm)"]
splitting_txts_ih = ["?", "#Deltam_{31}^{2} (atm)", "#Deltam_{21}^{2} (sol)"]

ys_nh.sort()
ys_ih.sort()

for i in xrange(3):

    dx_nh = (2.0/48, 4.0/48, 3.0/48)[i]
    dx_ih = (i+2.0)/48

    if i == 0:
        y0_nh = ax_y
        y0_ih = ax_y
    else:
        y0_nh = ys_nh[i-1]
        y0_ih = ys_ih[i-1]

    splitting_nh = ROOT.TPaveText(x_nh + dx_nh, y0_nh + 2*line_w,
                                  5.0/12, ys_nh[i] - 2*line_w)
    splitting_ih = ROOT.TPaveText(x_ih + dx_ih, y0_ih + 2*line_w,
                                  10.0/12, ys_ih[i] - 2*line_w)
    
    splitting_nh.AddText(splitting_txts_nh[i])
    splitting_ih.AddText(splitting_txts_ih[i])

    splittings.append(splitting_nh)
    splittings.append(splitting_ih)
    

for splitting in splittings:
    splitting.SetFillColor(ROOT.kWhite)
    splitting.SetBorderSize(0)
    splitting.SetTextSize(0.05)
    splitting.SetTextAlign(12)
    splitting.Draw()

# o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-
# Now draw the axes

ax_x = 3.0/48
l_zero = ROOT.TLine(ax_x, ax_y, 1-ax_x, ax_y)
l_zero.Draw()

ar_ax_x = 1.0/12
ar_size = 0.015
ar_l = ROOT.TArrow(ar_ax_x, 0, ar_ax_x, 0.9, ar_size, "|>")
ar_l.Draw()
ar_r = ROOT.TArrow(1-ar_ax_x, 0, 1-ar_ax_x, 0.9, ar_size, "|>")
ar_r.Draw()

# o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-o_o-
# Tick marks and labels for the axes

ticks = []
labels = []
for i in xrange(3):
    dx = 1.0/60
    dy = 0.05

    tick_nh = ROOT.TLine(ar_ax_x-dx, ys_nh[i],
                         ar_ax_x+dx, ys_nh[i])
    tick_ih = ROOT.TLine(1-(ar_ax_x-dx), ys_ih[i],
                         1-(ar_ax_x+dx), ys_ih[i])
    
    ticks.append(tick_nh)
    ticks.append(tick_ih)

    label_nh = ROOT.TPaveText(0, ys_nh[i] - dy,
                              ar_ax_x - dx, ys_nh[i] + dy)
    label_ih = ROOT.TPaveText(1-(ar_ax_x-dx), ys_ih[i] - dy,
                              1, ys_ih[i] + dy)

    label_nh.AddText("m_{%i}^{2}"%(i+1))
    label_ih.AddText("m_{%i}^{2}"%((3,1,2)[i]))

    labels.append(label_nh)
    labels.append(label_ih)


label_l = ROOT.TPaveText(0, ax_y - dy,
                         ar_ax_x - dx, ax_y + dy)
label_r = ROOT.TPaveText(1-(ar_ax_x-dx), ax_y - dy,
                         1, ax_y + dy)
label_l.AddText("0")
label_r.AddText("0")
labels.append(label_l)
labels.append(label_r)

label_topl = ROOT.TPaveText(ar_ax_x - dx, 0.9,
                            ar_ax_x + dx, 1)
label_topr = ROOT.TPaveText(1-(ar_ax_x+dx), 0.9,
                            1-(ar_ax_x-dx), 1)
label_topl.AddText("m^{2}")
label_topr.AddText("m^{2}")
labels.append(label_topl)
labels.append(label_topr)

label_normal = ROOT.TPaveText(2.0/12, 0.8,
                              5.0/12, 0.9)
label_inverted = ROOT.TPaveText(7.0/12, 0.8,
                                10.0/12, 0.9)
label_normal.AddText("\"normal\"")
label_normal.AddText("hierarchy")
label_inverted.AddText("\"inverted\"")
label_inverted.AddText("hierarchy")
labels.append(label_normal)
labels.append(label_inverted)

for tick in ticks:
    tick.Draw()

for label in labels:
    label.SetFillColor(ROOT.kWhite)
    label.SetBorderSize(0)
    label.SetTextSize(0.05)
    label.SetTextAlign(22)
    label.Draw()

dummy = raw_input("Press Enter...")
