# -*- coding: utf-8 -*-
from math import sqrt, cos, sin, pi
from cmath import exp
import numpy as np
from planegeometry.geometry import Triangle, Inversion # https://github.com/stla/PyPlaneGeometry
import matplotlib.pyplot as plt
from pymobiushyperbolic import gyrosegment, gyrocentroid, gyromidpoint


def reflection(A, B, M):
    circ = Triangle(A, B, gyromidpoint(A,B)).circumcircle()
    iota = Inversion.on_circle(circ)
    return iota.invert(M)

n = 3
p = 7 # 1/n + 1/p < 1/2 !
d = sqrt(
  cos(pi/n + pi/p) * cos(pi/p)
  / (sin(2*pi/p) * sin(pi/n) + cos(pi/n + pi/p) * cos(pi/p))
)
Sommets = [d] + [d*exp(1j*2*k*pi/n) for k in range(n)[1:]] 
Sommets = np.array([[z.real, z.imag] for z in Sommets])

def draw_triangle(axes, triangle, color, n=50):
    path = (
        gyrosegment(triangle.A, triangle.B, n=n)[:(n-1)]
        + gyrosegment(triangle.B, triangle.C, n=n)[:(n-1)]
        + gyrosegment(triangle.C, triangle.A, n=n)[:(n-1)]
    )
    axes.add_artist(
        plt.Polygon(
            path, closed=True, fill=True, facecolor=color, linewidth=2
        )
    )

figure, axes = plt.subplots(figsize=(10, 10))
axes.set_aspect(1)

def pavage(triangle, symetrie, niveau, i=1):
    global Centroids
    m = len(Centroids)
    rounded_centroids = np.vectorize(lambda x: round(x, 6), otypes=[float])(Centroids)
    _, indices = np.unique(rounded_centroids, return_index=True, axis=0)
    if len(indices) < m:
        Centroids = Centroids[indices, ]
    else:
        color = "navy" if i==1 else "yellow"
        draw_triangle(axes, triangle, color)
    if niveau > 0:
        for k in range(n):
            if k != symetrie:
                kp1 = (k+1) % n
                R = lambda M: reflection(Sommets[k], Sommets[kp1], M)
                newtriangle = Triangle(R(triangle.A), R(triangle.B), R(triangle.C))
                Centroids = np.vstack(
                    (
                        Centroids,
                        gyrocentroid(newtriangle.A, newtriangle.B, newtriangle.C)
                    )    
                )
                pavage(newtriangle, k, niveau-1, -i)

O = (0,0)
depth = 9
Centroids = np.empty((0,2), dtype=float)
for i in range(n):
    ip1 = (i+1) % n
    pavage(Triangle(O, Sommets[i], gyromidpoint(Sommets[i],Sommets[ip1])), 0, depth)
    im1 = n-1 if i == 0 else i-1
    pavage(Triangle(O, Sommets[i], gyromidpoint(Sommets[i],Sommets[im1])), 0, depth, -1)
plt.xlim(-1, 1)
plt.ylim(-1, 1)
plt.axis("off")
plt.show()
