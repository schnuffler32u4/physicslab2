import numpy as np
from functions import M
from scipy.optimize import curve_fit, least_squares
import matplotlib.pyplot as plt
import functions

V0 = 10.3
xdata, ydata = np.loadtxt("fittry.csv", delimiter=",", skiprows=1, unpack=True, encoding='UTF8')
xdata = 1000 * xdata * 2 * np.pi
ydata = ydata * np.sqrt(2)


# M = 0.06


def func(x, a, b, c, d, e):
    s = a * x ** 4 - M ** 2 * x ** 4 - c * x ** 2 + e
    t = d * x - b * x ** 3
    p = M * x ** 2 * V0
    A = p * t / (s ** 2 + t ** 2)
    B = p * s / (s ** 2 + t ** 2)
    return 1e10 * np.sqrt(A ** 2 + B ** 2)


def funcvar(x, L1, L2, R1, R2, C1):
    a = L1 * L2
    b = L1 * R2 + R1 * L2
    c = L1 / 1e-10 + R1 * R2 + L2 / C1
    d = R1 / 1e-10 + R2 / C1
    e = 1 / (C1 * 1e-10)
    return func(x, a, b, c, d, e)


popt, pcov = curve_fit(func, xdata, ydata, p0=[0.01, 1.4, 2e9, 1.4e14, 1e20], maxfev=500000)
print(M)

print(popt)
print(functions.circuit_values(*popt))

plt.plot(xdata, func(xdata, *popt), label='fit')
plt.plot(xdata, ydata)
plt.plot(xdata, funcvar(xdata, *functions.circuit_values(*popt)), label='yummy')
plt.legend()
plt.show()
