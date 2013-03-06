#!/bin/env python

from math import sin, cos, pow, exp, log, pi
import scipy.integrate

def ProjectedArea(theta, phi):
    # projected area for a particle incident from theta, phi
    # in cm^2
    r = 22.875-0137
    h = 2*20.44065
    return ((pi*r**2*fabs(cos(phi))*sin(theta) +
             2*r*h*sqrt(1-sin(theta)**2*cos(phi)**2)))

def flux(theta, h):

    log_intensity = (-(8e-4*h)/cos(theta)+(8e-4*h))
    
    intensity = exp(log_intensity)
    intensity *= pow(cos(theta),1.53)
    return intensity

def flux2d(theta, phi, h):
    return flux_sin(theta, h)

def vertical_flux(theta, h):
    return flux(theta, h)*cos(theta)

def bad_vertical_flux(theta, h):
    return flux(theta, h)*cos(theta)**2*sin(theta)

def vertical_flux2d(theta, phi, h):
    return vertical_flux_sin(theta, h)

def flux_sin(theta, h):
    intensity = flux(theta, h)*sin(theta)
    return intensity

def vertical_flux_sin(theta, h):
    intensity = vertical_flux(theta, h)*sin(theta)
    return intensity

def flux_sincos(theta, h):
    intensity = flux_sin(theta, h)*cos(theta)
    return intensity

def vertical_flux_sincos(theta, h):
    intensity = vertical_flux_sin(theta, h)*cos(theta)
    return intensity

def upper_list(pts):
    spts = sorted(pts)
    list = []

    # grab the upper most point on the left boundary
    mid_pt = spts[0]
    for pt in spts:
        if pt[0] == mid_pt[0]:
            mid_pt = pt
        else:
            break
        
    mid_y = mid_pt[1]
    for (x, y) in spts:
        if y >= mid_y:
            list.append((x,y))
    return list

def lower_list(pts):
    spts = sorted(pts)
    list = []

    # grab the lower most point on the left boundary
    mid_pt = spts[0]

    mid_y = spts[0][1]
    for (x, y) in spts:
        if y <= mid_y:
            list.append((x,y))
    return list

def boundary(x, pts):
    n = len(pts)
    if (x <= pts[0][0]):
        return pts[0][1]
    for i in xrange(n-1):
        x0 = pts[i][0]
        x1 = pts[i+1][0]
        if x0 < x <= x1:
            l = (x - x0)/(x1 - x0)
            y0 =  pts[i][1]
            y1 =  pts[i+1][1]
            y = y0 + l*(y1-y0)
            if y > pi/3:
                #return pi/3
                return y
            else:
                return y
    return pts[n-1][1]


good_phi = [0, 0.15, 0.9, 1.2, 1.2, 0.8]
good_phi.extend([-x for x in reversed(good_phi[1:])])
#good_phi.append(good_phi[0])
good_theta = [0.4, 0.45, 0.45, 0.75, 0.95, pi/3]
good_theta.extend(reversed(good_theta[1:]))
#good_theta.append(good_theta[0])

#good_phi.extend([-x for x in reversed(good_phi[:-1])])
#good_theta.extend(reversed(good_theta[:-1]))

good_pts = zip(good_phi, good_theta)
upper_boundary_pts = upper_list(good_pts)
lower_boundary_pts = lower_list(good_pts)

print(good_pts)
print(upper_boundary_pts)
print(lower_boundary_pts)

def upper_boundary(x):
    global upper_boundary_pts
    return boundary(x, upper_boundary_pts)

def lower_boundary(x):
    global lower_boundary_pts
    return boundary(x, lower_boundary_pts)

low_phi = min(upper_boundary_pts[0][0],
              lower_boundary_pts[0][0])
hi_phi = max(upper_boundary_pts[-1][0],
             lower_boundary_pts[-1][0])
print("integrating phi from %f to %f"%(low_phi, hi_phi))
            
base_depth = 1585
for depth in [x*base_depth for x in [0.9, 0.95, 1, 1.05, 1.1]]:

    print(" --- depth = %.3f --- "%(depth))

    i = (scipy.integrate.quad(bad_vertical_flux, 0, pi/2, (depth))[0]
         /scipy.integrate.quad(flux_sincos, 0, pi/2, (depth))[0])
    #i = (scipy.integrate.quad(bad_vertical_flux, 0, pi/2, (depth))[0])
    i_sin = (scipy.integrate.quad(vertical_flux_sin, 0, pi/2, (depth))[0]
             /scipy.integrate.quad(flux_sin, 0, pi/2, (depth))[0])
## i_sincos = (scipy.integrate.quad(vertical_flux_sincos, 0, pi/2, (depth))[0]
##             /scipy.integrate.quad(flux_sincos, 0, pi/2, (depth))[0])

## print("For 0 to pi/2:")
    print("  Esch ratio with cos?         : %0.2e"%(i))
    print("  Esch ratio without cos?      : %0.2e"%(i_sin))
## print("  Flux ratio (2 pi cos sin): %0.2e"%(1/i_sincos))

## i = (scipy.integrate.quad(vertical_flux, 0, pi/3, (depth))[0]
##      /scipy.integrate.quad(flux, 0, pi/3, (depth))[0])
## i_sin = (scipy.integrate.quad(vertical_flux_sin, 0, pi/3, (depth))[0]
##          /scipy.integrate.quad(flux_sin, 0, pi/3, (depth))[0])
## i_sincos = (scipy.integrate.quad(vertical_flux_sincos, 0, pi/3, (depth))[0]
##          /scipy.integrate.quad(flux_sincos, 0, pi/3, (depth))[0])

    I_h0 = 1/(2*pi*scipy.integrate.quad(flux_sin, 0, pi/2, (depth))[0])
    print("I(h, 0) = (Total Flux) * %0.2e"%(I_h0))

    I_h0_v = 1/(2*pi*scipy.integrate.quad(flux_sincos, 0, pi/2, (depth))[0])
    print("I(h, 0) = (Vertical Flux) * %0.2e"%(I_h0_v))
    
    I_h0_region = 1/(2*scipy.integrate.dblquad(flux2d,
                                               low_phi,
                                               hi_phi,
                                               lower_boundary,
                                               upper_boundary,
                                               (depth,))[0])
    print("I(h, 0) = (Flux in region) * %0.2e"%(I_h0_region))

## print("For 0 to pi/3:")
## print("  Flux ratio (nothing)     : %0.2e"%(1/i))
## print("  Flux ratio (2 pi sin)    : %0.2e"%(1/i_sin))
## print("  Flux ratio (2 pi cos sin): %0.2e"%(1/i_sincos))

    sa_region = 2*scipy.integrate.dblquad(lambda x,y: sin(x),
                                          low_phi,
                                          hi_phi,
                                          lower_boundary,
                                          upper_boundary)[0]
    sa_region = 1

    ratio_region = (2*scipy.integrate.dblquad(flux2d,
                                              low_phi,
                                              hi_phi,  
                                              lower_boundary,
                                              upper_boundary,
                                              (depth,))[0]
                    /(2*pi*scipy.integrate.quad(flux_sin, 0, pi/2, (depth))[0]))
    print("(Flux in region)/(Total flux up to pi/2): %0.3e"%(ratio_region/sa_region))

## ratio_region = (2*scipy.integrate.dblquad(flux2d,
##                                         good_pts[0][0],
##                                         good_pts[-1][0],
##                                         lower_boundary,
##                                         upper_boundary,
##                                         (depth,))[0]
##                 /(2*pi*scipy.integrate.quad(flux_sin, 0, pi/3, (depth))[0]))
## print("(Flux in region)/(Total flux up to pi/3): %0.2e"%(ratio_region))

## ratio_region = (2*scipy.integrate.dblquad(vertical_flux2d,
##                                         good_pts[0][0],
##                                         good_pts[-1][0],
##                                         lower_boundary,
##                                         upper_boundary,
##                                         (depth,))[0]
##                 /(2*pi*scipy.integrate.quad(vertical_flux_sin,
##                                             0, pi/2, (depth))[0]))
## print("(Vertical flux in region)/(Vertical flux up to pi/2): %0.2e"%(ratio_region))

## ratio_region = (2*scipy.integrate.dblquad(vertical_flux2d,
##                                         good_pts[0][0],
##                                         good_pts[-1][0],
##                                         lower_boundary,
##                                         upper_boundary,
##                                         (depth,))[0]
##                 /(2*pi*scipy.integrate.quad(vertical_flux_sin,
##                                             0, pi/3, (depth))[0]))
## print("(Vertical flux in region)/(Vertical flux up to pi/3): %0.2e"%(ratio_region))

    ratio_region = (2*scipy.integrate.dblquad(flux2d,
                                              low_phi,
                                              hi_phi,
                                              lower_boundary,
                                              upper_boundary,
                                              (depth,))[0]
                    /(2*pi*scipy.integrate.quad(vertical_flux_sin,
                                                0, pi/2, (depth))[0]))
    print("(Flux in region)/(Vertical flux up to pi/2): %0.3e"%(ratio_region/sa_region))

## ratio_region = (2*scipy.integrate.dblquad(flux2d,
##                                         good_pts[0][0],
##                                         good_pts[-1][0],
##                                         lower_boundary,
##                                         upper_boundary,
##                                         (depth,))[0]
##                 /(2*pi*scipy.integrate.quad(vertical_flux_sin,
##                                             0, pi/3, (depth))[0]))
## print("(Flux in region)/(Vertical flux up to pi/3): %0.2e"%(ratio_region))
    
