# -*- coding: utf-8 -*-
from math import sqrt
import numpy as np
import quaternion
import pyvista as pv

# plane passing by points p1, p2, p3 #
def plane3pts(p1, p2, p3):
    xcoef = (p1[1] - p2[1]) * (p2[2] - p3[2]) - (p1[2] - p2[2]) * (p2[1] - p3[1])
    ycoef = (p1[2] - p2[2]) * (p2[0] - p3[0]) - (p1[0] - p2[0]) * (p2[2] - p3[2])
    zcoef = (p1[0] - p2[0]) * (p2[1] - p3[1]) - (p1[1] - p2[1]) * (p2[0] - p3[0])
    offset = p1[0] * xcoef + p1[1] * ycoef + p1[2] * zcoef
    return np.array([xcoef, ycoef, zcoef, offset])


# center, radius and normal of the circle passing by three points #
def circleCenterAndRadius(p1, p2, p3):
    p12 = (p1 + p2) / 2
    p23 = (p2 + p3) / 2
    v12 = p2 - p1
    v23 = p3 - p2
    plane = plane3pts(p1, p2, p3)
    A = np.column_stack((plane[0:3], v12, v23))
    b = np.array([plane[3], np.vdot(p12, v12), np.vdot(p23, v23)])
    center = np.matmul(np.linalg.inv(np.transpose(A)), b)
    r = np.linalg.norm(p1 - center)
    return dict(center=center, radius=r, normal=plane[0:3])


# map half-space to unit ball (Fenchel p.52) ####
# rk: x3 = 0 <=> Mod(Phi(x))=1
def Phi(x1x2x3):
    q = quaternion.from_float_array(x1x2x3 + [0])
    Hj = quaternion.y
    img = -Hj - 2*(q+Hj).inverse()
    return np.array([img.w, img.x, img.y])

def PhiInv(x1x2x3):
    q = quaternion.from_float_array(x1x2x3 + [0])
    Hj = quaternion.y
    img = (-0.5*(q + Hj)).inverse() - Hj
    return np.array([img.w, img.x, img.y])

# draw hyperbolic polyhedron ####
# helper function
def f(p1, p2):
    return [(p1[0]+p2[0])/2, (p1[1]+p2[1])/2, np.linalg.norm(p1-p2)/2]

def hpolyhedron(pltr, vertices, edges, col = "black"):
    P = [PhiInv(vertex) for vertex in vertices]
    for edge in edges: 
        i1 = edge[0]; i2 = edge[1]
        Q = Phi(f(P[i1], P[i2]))
        vertices = np.array(vertices)
        Q1 = vertices[i1, :]; Q2 = vertices[i2, :]
        circle = circleCenterAndRadius(Q1, Q2, Q)
        arc = pv.CircularArc(Q1, Q2, circle["center"])
        pltr.add_mesh(arc.tube(radius=0.02), color=col, smooth_shading=True)
    pltr.add_points(vertices, render_points_as_spheres=True, point_size=10)

# vertices ####
phi = (1+sqrt(5))/2
a = 1/sqrt(3)
b = a/phi 
c = a*phi
vertices = [
    [ a,  a,  a], 
    [ a,  a, -a],
    [ a, -a,  a],
    [-a, -a,  a],
    [-a,  a, -a],
    [-a,  a,  a],
    [ 0,  b, -c], 
    [ 0, -b, -c], 
    [ 0, -b,  c],
    [ c,  0, -b],
    [-c,  0, -b],
    [-c,  0,  b],
    [ b,  c,  0],
    [ b, -c,  0],
    [-b, -c,  0],
    [-b,  c,  0],
    [ 0,  b,  c],
    [ a, -a, -a],
    [ c,  0,  b],
    [-a, -a, -a]
  ]

edges = [
 [0, 12],
 [0, 16],
 [0, 18],
 [1, 6],
 [1, 9],
 [1, 12],
 [2, 8],
 [2, 13],
 [2, 18],
 [3, 8],
 [3, 11],
 [3, 14],
 [4, 6],
 [4, 10],
 [4, 15],
 [5, 11],
 [5, 15],
 [5, 16],
 [6, 7],
 [7, 17],
 [7, 19],
 [8, 16],
 [9, 17],
 [9, 18],
 [10, 11],
 [10, 19],
 [12, 15],
 [13, 14],
 [13, 17],
 [14, 19]
]


pltr = pv.Plotter()
hpolyhedron(pltr, vertices, edges, col="navy")
pltr.show()
