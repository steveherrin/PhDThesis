#!/bin/env python

import random
import numpy
from math import log, exp, cos, sin, pi

panel_width = 76.2
panel_length = 305
panel_depth = 2.54
panel_sep = 30.5-1*panel_depth

def intensity(theta, depth):
    log_i = (-(8e-4*depth)/cos(theta)+
                     1.53*log(cos(theta))+(8e-4*depth))
    i = sin(theta)*exp(log_i)
    #i = exp(log_i)
    #i = sin(theta)*pow(cos(theta),1.53)*exp(-8e-4*depth*(1/cos(theta)-1))
    return i

def monte_carlo(N, depth):
    global panel_width
    global panel_length
    global panel_depth
    global panel_dep

    hit = 0
    max_intensity = 0
    seen_at = 0
    for theta in numpy.linspace(0, pi/2, 1e6, False):
        i = cos(theta)*intensity(theta, depth)
        if i > max_intensity:
            max_intensity = i
            at = theta

    print("Maximum intensity for %f m.w.e. is %0.4e"%(depth, max_intensity))

    for i in xrange(N):

        while(True):
            theta = pi/2*random.random()

            y = max_intensity*random.random()

            int = cos(theta)*intensity(theta, depth)
            if (int >= y):
                break

        phi = 2*pi*(random.random()-0.5)

        x0 = random.random()*panel_width
        y0 = random.random()*panel_length
        z0 = panel_sep + 2*panel_depth

        # z0 = 0 so
        r = z0/cos(theta)
        x1 = x0 - r*cos(phi)*sin(theta)
        y1 = y0 - r*sin(phi)*sin(theta)

        # also check the upper surface of the lower panel
        r = (z0-panel_depth)/cos(theta)
        x2 = x0 - r*cos(phi)*sin(theta)
        y2 = y0 - r*sin(phi)*sin(theta)

        if (((0 <= x1 < panel_width) and (0 <= y1 < panel_length))
            or ((0 <= x2 < panel_width) and (0 <= y2 < panel_length))):
            hit += 1
            
    print("For depth %0.1f m.w.e., %i hits out of %i tries (= %0.1f %%)"
          %(depth, hit, N, 100.0*hit/N))

    return hit/N


#monte_carlo(1000000, 1585)

for depth in [x*1585 for x in [0.9, 1.0, 1.1]]:
    monte_carlo(1000000, depth)
