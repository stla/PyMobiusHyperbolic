__version__ = '0.1.0'

from math import sqrt, tanh, atanh
import numpy as np
import pyvista as pv


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
    return gyroscalar(
        0.5, (gA2*A + gB2*B + gC2*C) / (gA2 + gB2 + gC2 - 1.5), s=s
    )

def gyrotube(A, B, s, r, npoints=300):
    """
    Tubular hyperbolic segment.

    Parameters
    ----------
    A,B : points (lists or arrays)
        The two endpoints of the segment.
    s : positive float
        Curvature parameter.
    r : positive float
        Radius of the tube.
    npoints : integer
        Number of points along the segment. The default is 300.

    Returns
    -------
    PyVista mesh
        A PyVista mesh ready for inclusion in a plotting region.

    """
    AB = gyrosegment(A, B, n=100, s=s)
    splineAB = pv.Spline(AB, npoints)
    splineAB["scalars"] = np.arange(splineAB.n_points)
    return splineAB.tube(radius=r).extract_geometry()


def gyrosubdiv(A1, A2, A3, s):
    M12 = gyromidpoint(A1, A2, s)
    M13 = gyromidpoint(A1, A3, s)
    M23 = gyromidpoint(A2, A3, s)
    return [[A1, M12, M13], [A2, M23, M12], [A3, M13, M23], [M12, M13, M23]]


def gyrotriangle(A, B, C, s, depth=5, tol=1e-6):
    """
    Hyperbolic triangle.

    Parameters
    ----------
    A,B,C : points (lists or arrays)
        The vertices of the triangle.
    s : positive float
        Curvature parameter.
    depth : integer
        The number of recursive subdivions. The default is 5.
    tol : small positive float
        The tolerance used to merge duplicated points in the mesh.
        The default is 1e-6.

    Returns
    -------
    PyVista mesh
        A PyVista mesh ready for inclusion in a plotting region.

    """
    subd = gyrosubdiv(A, B, C, s)
    for _ in range(depth - 1):
        lst = list(map(lambda t: gyrosubdiv(*t, s), subd))
        subd = np.concatenate(lst)
    tsubd = list(map(np.transpose, subd))
    vertices = np.transpose(np.concatenate(tsubd, axis=1))
    nvertices = vertices.shape[0]
    ntriangles = nvertices // 3
    repeats3 = np.full((ntriangles, 1), 3)
    indices = np.array(range(nvertices)).reshape(ntriangles, 3)
    faces = np.hstack((repeats3, indices)).flatten()
    return pv.PolyData(vertices, faces).clean(tolerance=tol)
