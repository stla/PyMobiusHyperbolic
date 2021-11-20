# -*- coding: utf-8 -*-
from math import sqrt
import numpy as np
import pyvista as pv
from pymobiushyperbolic import gyrotube, gyrotriangle

####~~ bilunabirotunda ~~####

def hpolyhedron(pltr, vertices, faces, edges, s):
    for face in faces:
        vs = vertices[face, :]
        mesh = gyrotriangle(vs[0], vs[1], vs[2], depth=4, s=s)
        pltr.add_mesh(mesh, smooth_shading=True, color="navy")
    for edge in edges:
        A = vertices[edge[0], :]
        B = vertices[edge[1], :]
        mesh = gyrotube(A, B, s=s, r=0.035)
        pltr.add_mesh(mesh, smooth_shading=True, color="yellow")
    for vertex in vertices:
        sphere = pv.Sphere(0.05, center=vertex)
        pltr.add_mesh(sphere, smooth_shading=True, color="yellow")

# vertices ####
phi = (1+sqrt(5))/2
vertices = np.array([ 
    [1, 0, phi**2],
    [-1, 0, phi**2],
    [1, 0, -phi**2],
    [-1, 0, -phi**2],
    [phi, 1, 1],
    [-phi, 1, 1],
    [phi, -1, 1],
    [phi, 1, -1],
    [-phi, -1, 1],
    [-phi, 1, -1],
    [phi, -1, -1],
    [-phi, -1, -1],
    [0, phi ,0],
    [0, -phi, 0]
  ])
M = np.max(np.linalg.norm(vertices, axis=1))
s = M * 1.02

faces = [
 [7, 12, 4],
 [7, 10, 2],
 [9, 12, 5],
 [6, 10, 13],
 [6, 0, 4],
 [11, 9, 3],
 [8, 1, 5],
 [8, 11, 13],
 [1, 0, 4],
 [1, 12, 4],
 [1, 12, 5],
 [9, 3, 2],
 [9, 7, 12],
 [9, 7, 2],
 [6, 7, 10],
 [6, 7, 4],
 [11, 10, 2],
 [11, 3, 2],
 [11, 10, 13],
 [8, 11, 9],
 [8, 9, 5],
 [8, 6, 0],
 [8, 1, 0],
 [8, 6, 13]
]

edges = [
 [0, 1],
 [0, 4],
 [0, 6],
 [1, 5],
 [1, 8],
 [2, 3],
 [2, 7],
 [2, 10],
 [3, 9],
 [3, 11],
 [4, 6],
 [4, 7],
 [4, 12],
 [5, 8],
 [5, 9],
 [5, 12],
 [6, 10],
 [6, 13],
 [7, 10],
 [7, 12],
 [8, 11],
 [8, 13],
 [9, 11],
 [9, 12],
 [10, 13],
 [11, 13]
]

pltr = pv.Plotter()
hpolyhedron(pltr, vertices, faces, edges, s)
pltr.show()
