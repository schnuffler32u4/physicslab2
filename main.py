import numpy as np
from functions import M
from scipy.optimize import curve_fit, least_squares
import matplotlib.pyplot as plt

V0 = 10.3
xdata, ydata = np.loadtxt("fittry.csv", delimiter=",", skiprows=1, unpack=True, encoding='UTF8')
xdata = 1000 * xdata / (2 * np.pi)
ydata = ydata * np.sqrt(2)
M = 0.09


def func(x, a, b, c, d, e):
    s = a * x ** 4 - M ** 2 * x ** 4 - c * x ** 2 + e
    t = d * x - b * x ** 3
    p = M * x ** 2 * V0
    A = p * t / (s ** 2 + t ** 2)
    B = p * s / (s ** 2 + t ** 2)
    return 1e10 * np.sqrt(A ** 2 + B ** 2)


popt, pcov = curve_fit(func, xdata, ydata, p0=[0.01, 1.4, 2e9, 1.4e14, 1e20],  maxfev=500000)
print(M)
# x = np.linspace(0, 900000, 10000)
print(popt)
# print(pcov)
plt.plot(xdata, func(xdata, *popt), label='fit')
plt.plot(xdata, ydata)
# plt.plot(x, func(x, 0.01, 1.4, 2e9, 1.4e11, 1e20))
plt.legend()
plt.show()
