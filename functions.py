import scipy.integrate as si
import numpy as np
from scipy.optimize import fsolve


def mutual(dis):
    """The function takes the distance in meters as an input and outputs the mutual inductance between the coils."""
    R = 0.105

    def f(x, y):
        return 200 ** 2 * 10 ** (-7) * R ** 2 * np.cos(x - y) / (np.sqrt(2 * R ** 2 * (1 - np.cos(x - y)) + dis ** 2))

    return si.dblquad(f, 0, 2 * np.pi, 0, 2 * np.pi)[0]


def circuit_values(a, b, c, d, e):
    """The function takes as input the values of a,b,c,d,e and returns the corresponding capacitances, resistances
    and so on. We assume C2 = 1e-10. The output is an array of the form [L1,L2,R1,R2,C1]"""
    C2 = 1e-10
    C1 = 1 / (e * C2)

    def leq(x):
        return a / (x * C2) + C1 * (b - a * C1 * d /x ) / (x - a * C1 / (x * C2)) * (
                d - (b - a * C1 * d / x) / (C2 * (x - a * C1 / (x * C2)))) + x / C1 - c

    L2 = fsolve(leq, 1, maxfev=500000)[0]
    L1 = a / L2
    R1 = (b - a * C1 * d / L2) / (L2 - a * C1 / (L2 * C2))
    R2 = C1 * (d - R1 / C2)
    return [L1, L2, R1, R2, C1]
