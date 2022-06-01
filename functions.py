import scipy.integrate as si
import numpy as np


def midpoint_double(f, a, b, c, d, nx, ny):
    hx = (b - a) / nx
    hy = (d - c) / ny
    I = 0
    for i in range(nx):
        for j in range(ny):
            xi = a + hx / 2 + i * hx
            yj = c + hy / 2 + j * hy
            I += hx * hy * f(xi, yj, dis)
    return I


R = 0.105
dis = 0.05


def f(x, y):
    return 200 ** 2 * 10 ** (-7) * R ** 2 * np.cos(x - y) / (np.sqrt(2 * R ** 2 * (1 - np.cos(x - y)) + dis ** 2))

M = si.dblquad(f, 0, 2*np.pi, 0, 2*np.pi)[0]



