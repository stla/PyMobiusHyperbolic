# -*- coding: utf-8 -*-
from math import sqrt
import numpy as np
import pyvista as pv
from pymobiushyperbolic import gyrotube, gyrotriangle

####~~ great stellated dodecahedron ~~####

def hpolyhedron(pltr, vertices, faces, edges, s):
    for face in faces:
        vs = vertices[face, :]
        mesh = gyrotriangle(vs[0], vs[1], vs[2], depth=5, s=s)
        pltr.add_mesh(mesh, smooth_shading=True, color="navy")
    for edge in edges:
        A = vertices[edge[0], :]
        B = vertices[edge[1], :]
        mesh = gyrotube(A, B, s=s, r=0.015)
        pltr.add_mesh(mesh, smooth_shading=True, color="yellow")
    for vertex in vertices:
        sphere = pv.Sphere(0.025, center=vertex)
        pltr.add_mesh(sphere, smooth_shading=True, color="yellow")

# vertices ####
vertices = np.array([ 
    [0, 0, (1+sqrt(5))/sqrt(3)*3/2], 
    [(1+sqrt(5))/sqrt(3), 0, -(1+sqrt(5))*sqrt(5/3)/2], 
    [-(1+sqrt(5))/sqrt(3)/2, (1+sqrt(5))/2, -(1+sqrt(5))*sqrt(5/3)/2], 
    [-(1+sqrt(5))/sqrt(3)/2, -(1+sqrt(5))/2, -(1+sqrt(5))*sqrt(5/3)/2],
    [-(1+sqrt(5))*sqrt(5/3)/2, (1+sqrt(5))/2, (1+sqrt(5))/sqrt(3)/2],
    [-(1+sqrt(5))*sqrt(5/3)/2, -(1+sqrt(5))/2, (1+sqrt(5))/sqrt(3)/2],
    [-1/(1+sqrt(5))/2/sqrt(3), -(sqrt(5)+3)/2, (1+sqrt(5))/sqrt(3)/2], 
    [(2+sqrt(5))/sqrt(3), -1, (1+sqrt(5))/sqrt(3)/2], 
    [(2+sqrt(5))/sqrt(3), 1, (1+sqrt(5))/sqrt(3)/2],
    [-1/(1+sqrt(5))/2/sqrt(3), (sqrt(5)+3)/2, (1+sqrt(5))/sqrt(3)/2],
    [1/(1+sqrt(5))/2/sqrt(3), -(sqrt(5)+3)/2, -(1+sqrt(5))/sqrt(3)/2],
    [1/(1+sqrt(5))/2/sqrt(3), (sqrt(5)+3)/2, -(1+sqrt(5))/sqrt(3)/2],
    [(1+sqrt(5))*sqrt(5/3)/2, (1+sqrt(5))/2, -(1+sqrt(5))/sqrt(3)/2],
    [-(2+sqrt(5))/sqrt(3), -1, -(1+sqrt(5))/sqrt(3)/2],
    [-(2+sqrt(5))/sqrt(3), 1, -(1+sqrt(5))/sqrt(3)/2],
    [(1+sqrt(5))*sqrt(5/3)/2, -(1+sqrt(5))/2, -(1+sqrt(5))/sqrt(3)/2],
    [(1+sqrt(5))/sqrt(3)/2, (1+sqrt(5))/2, (1+sqrt(5))*sqrt(5/3)/2],
    [(1+sqrt(5))/sqrt(3)/2, -(1+sqrt(5))/2, (1+sqrt(5))*sqrt(5/3)/2],
    [-(1+sqrt(5))/sqrt(3), 0, (1+sqrt(5))*sqrt(5/3)/2],
    [0, 0, -(1+sqrt(5))/sqrt(3)*3/2],
    [1/sqrt(15), 1/sqrt(5), (3/sqrt(5)-1)/sqrt(3)/2], 
    [-2/sqrt(15), 0, (3/sqrt(5)-1)/sqrt(3)/2], 
    [1/sqrt(15), -1/sqrt(5), (3/sqrt(5)-1)/sqrt(3)/2],
    [-(sqrt(5)-1)/sqrt(15), 0, -(1+sqrt(5))/2/sqrt(5)/sqrt(3)], 
    [(sqrt(5)-1)/sqrt(15)/2, -(5-sqrt(5))/10, -(1+sqrt(5))/2/sqrt(5)/sqrt(3)], 
    [(sqrt(5)-1)/sqrt(15)/2, (5-sqrt(5))/10, -(1+sqrt(5))/2/sqrt(5)/sqrt(3)],
    [-(sqrt(5)-1)/sqrt(15)/2, -(5-sqrt(5))/10, (1+sqrt(5))/2/sqrt(5)/sqrt(3)],
    [-(sqrt(5)-1)/sqrt(15)/2, (5-sqrt(5))/10, (1+sqrt(5))/2/sqrt(5)/sqrt(3)],
    [(sqrt(5)-1)/sqrt(15), 0, (1+sqrt(5))/2/sqrt(5)/sqrt(3)],
    [2/sqrt(15), 0, -(3/sqrt(5)-1)/sqrt(3)/2],
    [-1/sqrt(15), 1/sqrt(5), -(3/sqrt(5)-1)/sqrt(3)/2],
    [-1/sqrt(15), -1/sqrt(5), -(3/sqrt(5)-1)/sqrt(3)/2]
  ])
M = np.max(np.linalg.norm(vertices, axis=1))
vertices = vertices/M * 0.95


faces = [
 [20, 0, 1],
 [20, 1, 4],
 [20, 4, 7],
 [20, 7, 2],
 [20, 2, 0],
 [21, 0, 2],
 [21, 2, 6],
 [21, 6, 9],
 [21, 9, 3],
 [21, 3, 0],
 [22, 0, 3],
 [22, 3, 8],
 [22, 8, 5],
 [22, 5, 1],
 [22, 1, 0],
 [23, 1, 5],
 [23, 5, 11],
 [23, 11, 10],
 [23, 10, 4],
 [23, 4, 1],
 [24, 2, 7],
 [24, 7, 13],
 [24, 13, 12],
 [24, 12, 6],
 [24, 6, 2],
 [25, 3, 9],
 [25, 9, 15],
 [25, 15, 14],
 [25, 14, 8],
 [25, 8, 3],
 [26, 4, 10],
 [26, 10, 16],
 [26, 16, 13],
 [26, 13, 7],
 [26, 7, 4],
 [27, 5, 8],
 [27, 8, 14],
 [27, 14, 17],
 [27, 17, 11],
 [27, 11, 5],
 [28, 6, 12],
 [28, 12, 18],
 [28, 18, 15],
 [28, 15, 9],
 [28, 9, 6],
 [29, 10, 11],
 [29, 11, 17],
 [29, 17, 19],
 [29, 19, 16],
 [29, 16, 10],
 [30, 12, 13],
 [30, 13, 16],
 [30, 16, 19],
 [30, 19, 18],
 [30, 18, 12],
 [31, 14, 15],
 [31, 15, 18],
 [31, 18, 19],
 [31, 19, 17],
 [31, 17, 14]
]

edges = [
 [0, 16],
 [0, 17],
 [0, 18],
 [1, 12],
 [1, 15],
 [1, 19],
 [2, 11],
 [2, 12],
 [2, 14],
 [2, 19],
 [3, 10],
 [3, 13],
 [3, 15],
 [3, 19],
 [4, 9],
 [4, 11],
 [4, 14],
 [4, 16],
 [4, 18],
 [5, 6],
 [5, 10],
 [5, 13],
 [5, 17],
 [5, 18],
 [6, 7],
 [6, 10],
 [6, 15],
 [6, 17],
 [7, 8],
 [7, 15],
 [7, 17],
 [8, 9],
 [8, 12],
 [8, 16],
 [9, 11],
 [9, 12],
 [9, 16],
 [10, 13],
 [10, 15],
 [11, 12],
 [11, 14],
 [13, 14]
]

pltr = pv.Plotter()
hpolyhedron(pltr, vertices, faces, edges, s=1)
pltr.show()
