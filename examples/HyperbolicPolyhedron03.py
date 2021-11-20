# -*- coding: utf-8 -*-
from math import sqrt
import numpy as np
import pyvista as pv
from pymobiushyperbolic import gyrotube, gyrotriangle

####~~ icosahedron ~~####

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
 [0, 0.618033988749895, 1],
 [0, 0.618033988749895, -1],
 [0, -0.618033988749895, 1],
 [0, -0.618033988749895, -1],
 [0.618033988749895, 1, 0],
 [0.618033988749895, -1, 0],
 [-0.618033988749895, 1, 0],
 [-0.618033988749895, -1, 0],
 [1, 0, 0.618033988749895],
 [-1, 0, 0.618033988749895],
 [1, 0, -0.618033988749895],
 [-1, 0, -0.618033988749895]
])
# M = np.max(np.linalg.norm(vertices, axis=1))
# vertices = vertices/M * 0.95


faces = [
 [0, 2, 8],
 [0, 8, 4],
 [0, 4, 6],
 [0, 6, 9],
 [0, 9, 2],
 [3, 11, 1],
 [3, 1, 10],
 [3, 10, 5],
 [3, 5, 7],
 [3, 7, 11],
 [8, 2, 5],
 [4, 8, 10],
 [6, 4, 1],
 [9, 6, 11],
 [2, 9, 7],
 [1, 11, 6],
 [10, 1, 4],
 [5, 10, 8],
 [7, 5, 2],
 [11, 7, 9]
]

edges = [
 [0, 2],
 [1, 3],
 [0, 4],
 [1, 4],
 [2, 5],
 [3, 5],
 [0, 6],
 [1, 6],
 [4, 6],
 [2, 7],
 [3, 7],
 [5, 7],
 [0, 8],
 [2, 8],
 [4, 8],
 [5, 8],
 [0, 9],
 [2, 9],
 [6, 9],
 [7, 9],
 [1, 10],
 [3, 10],
 [4, 10],
 [5, 10],
 [8, 10],
 [1, 11],
 [3, 11],
 [6, 11],
 [7, 11],
 [9, 11]
]

pltr = pv.Plotter(window_size=[512,512], off_screen=False)
hpolyhedron(pltr, vertices/sqrt(2-2/(1+sqrt(5)))*0.95, faces, edges, s=1)
pltr.show()
