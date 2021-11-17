__version__ = '0.1.0'

from math import sqrt, tanh, atanh
import numpy as np


def dotprod(x, y=None):
    if y is None:
        y = x
    return np.vdot(x, y)   


def gyroadd(X, Y, s=1):
    s2 = s*s
    dX = dotprod(X)
    dY = dotprod(Y)
    dXY = dotprod(X, Y)
    return (
        (1 + 2*dXY/s2 + dY/s2) * X + (1 - dX/s2) * Y
    ) / (1 + 2*dXY/s2 + dX*dY/s2/s2)


def gyroscalar(r, X, s=1):
    Xnorm = sqrt(dotprod(X))
    return s/Xnorm * tanh(r * atanh(Xnorm/s)) * X


def gyrovector(A, B, s=1):
    return gyroadd(-A, B, s=s)


def gyroABt(A, B, t, s=1):
    return gyroadd(A, gyroscalar(t, gyrovector(A, B, s=s), s=s), s=s)


def gyrosegment(A, B, n=50, s=1):
    return [gyroABt(A, B, t, s=s) for t in np.linspace(0, 1, n)]


def gyromidpoint(A, B, s=1):
    return gyroABt(A, B, t=0.5, s=s)


def gyrocentroid(A, B, C, s=1):
    gA2 = 1/(1 - dotprod(A))
    gB2 = 1/(1 - dotprod(B))
    gC2 = 1/(1 - dotprod(C))
    return gyroscalar(0.5, (gA2*A + gB2*B + gC2*C) / (gA2 + gB2 + gC2 - 1.5), s=s)
