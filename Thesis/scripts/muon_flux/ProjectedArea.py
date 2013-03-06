#!/bin/env python

import numpy
from math import pi, sin, cos, fabs, sqrt
import sys
import time
import ROOT

def RotateCoordinates(points, theta, phi):
    # We want to be in a coordinate system such that
    # theta is the zenith angle
    # and phi is the angle from north (such that pi/2 is east)
    rot_phi = numpy.matrix([[ cos(phi), 0, -sin(phi)],
                            [ 0,        1,  0       ],
                            [ sin(phi), 0,  cos(phi)]])
    
    rot_theta = numpy.matrix([[1,  0,               0         ],
                              [0,  cos(pi/2-theta), sin(pi/2-theta)],
                              [0, -sin(pi/2-theta), cos(pi/2-theta)]])
    
    rotation = rot_theta*rot_phi
    
    rotated_points = []
    for point in points:
        v_point = numpy.matrix(point)
        if (numpy.shape(v_point)[0] != 3):
            v_point = v_point.transpose()
        rotated_point = rotation * v_point
        rotated_points.append(rotated_point)

    return rotated_points

def Test():
    phi = 0
    theta = pi/2
    print("Theta = %f, Phi = %f"%(theta, phi))

    points = [[0, 0, 1],
              [0, 0, -1],
              [0, 1, 0],
              [0, -1, 0],
              [1, 0, 0],
              [-1, 0, 0]]
    
    rotated_points = RotateCoordinates(points, theta, phi)
    for (point, rotated_point) in zip(points, rotated_points):
        print("---")
        print point
        print "  rotated to"
        print rotated_point


def IsRightTurn(a, b, c):
    # returns true if you can get from a to b to c with a right turn at b
    # and returns false if you must make a left turn

    # sign(sin(x)) = sign((b-a)x(c-a))
    # where x is the angle between b-a and c-a, and x is cross product
    sign = ((b[0]-a[0])*(c[1]-a[1])-(c[0]-a[0])*(b[1]-a[1]))
    return (sign < 0)


def KeepRightTurns(hull, new_point):
    # Try adding new point to the convex hull of the previous points
    # If the new point creates a left turn, then things aren't convex anymore,
    # so start chopping points off the current hull until it's all right turns
    while (len(hull) > 1 and not IsRightTurn(hull[-2], hull[-1], new_point)):
        hull.pop()
    if (len(hull) == 0) or (hull[-1] != new_point):
        hull.append(new_point)
    return hull


def GetConvexHull(points):
    # returns a list of points on the convex hull of a list of points
    # the points are returned in clockwise order
    points = sorted(points)
    lower_hull = []
    upper_hull = []
    for i in range(len(points)):
        lower_hull = KeepRightTurns(lower_hull, points[i])
        upper_hull = KeepRightTurns(upper_hull, points[-(i+1)])

    # now lower_hull connects bottom left to upper right
    # upper_hull connects upper right to bottom left
    # and they share the same endpoints, so link them up
    
    lower_hull.extend(upper_hull[i] for i in range(1, len(upper_hull)-1))
    return lower_hull


def GetArea(points):
    # get the area of the convex polygon from a list of points that are
    # the vertices of the polygon in order around it
    n = len(points)
    if n < 3:
        return 0
                      
    sum = points[n-1][0]*points[0][1] - points[0][0]*points[n-1][1]
    for i in range(n-1):
        sum += points[i][0]*points[i+1][1]
        sum -= points[i+1][0]*points[i][1]
    return 0.5*fabs(sum)


def CylToCartesian(r, z, phi):
    x = r*cos(phi)
    y = r*sin(phi)
    return [x, y, z]


def GetCylinderPoints():
    points = []
    max_z = 0.5
    max_r = 1
    for z in numpy.linspace(-max_z, max_z, 2, True):
        for phi in numpy.linspace(0, 2*pi, 10*360, False):
            points.append(CylToCartesian(max_r, z, phi))
    return points

def RotateAndProject(points, theta, phi):
    rotated = RotateCoordinates(points, theta, phi)
    projected = []
    for point in rotated:
        projected.append([point[0], point[1]])
    return projected

def ProjectedArea(theta, phi):
    points = RotateAndProject(GetCylinderPoints(), theta, phi)
    
    g_pts = ROOT.TGraph(len(points))
    for (i, point) in enumerate(points):
        g_pts.SetPoint(i, point[0], point[1])
    
    convex_hull = GetConvexHull(points)

    area = GetArea(convex_hull)
    
    g_hull = ROOT.TGraph(len(convex_hull)+1)
    for (i, point) in enumerate(convex_hull):
        g_hull.SetPoint(i, point[0], point[1])
    g_hull.SetPoint(len(convex_hull), convex_hull[0][0], convex_hull[0][1])
    g_hull.SetLineColor(ROOT.kRed)

    c = ROOT.TCanvas("c","Projection",480,480)
    g_pts.SetTitle("Projected Area = %0.2f"%(area))
    g_pts.Draw("AP")
    g_hull.Draw("L")

    #time.sleep(1)

    #print("Area is %f"%(area))
    #dummy = raw_input("Press Enter...")
    
    return area

if '__main__' in __name__:
    n_p = 48
    n_t = 24
    g_raw = ROOT.TGraph2D(n_p*n_t)
    g_raw.SetName("g_raw")
    g_sub = ROOT.TGraph2D(n_p*n_t)
    g_sub.SetName("g_sub")
    i = 0

    sub_file = ROOT.TFile("projected_areas_r1l2.root","read")
    g_old = sub_file.Get("g_raw")
    
    for theta in numpy.linspace(0, pi, n_t):
        for phi in numpy.linspace(0, 2*pi, n_p):
            area = ProjectedArea(theta, phi)
            g_raw.SetPoint(i, theta, phi, area)
            #g_sub.SetPoint(i, theta, phi, area - g_old.GetZ()[i])
            prediction = pi*fabs(cos(phi))*sin(theta)+2*sqrt(1-sin(theta)**2*fabs(cos(phi)**2))
            g_sub.SetPoint(i, theta, phi, area - prediction)
            i += 1
            if (i%10 == 0):
                print "Done with %i of %i"%(i,n_p*n_t)
    sub_file.Close()
                
    save_file = ROOT.TFile("projected_areas.root","recreate")
    g_raw.Write("g_raw")
    g_sub.Write("g_sub")
    save_file.Close()

    f = ROOT.TF2("f","TMath::Pi()*abs(cos(y))*sin(x)+2*sqrt(1-sq(sin(x))*sq(cos(y)))",0,pi,0,2*pi)

    save_file = ROOT.TFile("projected_areas.root","update")
    f.Write("f")
    save_file.Close()

    c1 = ROOT.TCanvas("c1")
    g_raw.Draw("CONT4Z")
    c2 = ROOT.TCanvas("c2")
    f.Draw("CONT4Z")
