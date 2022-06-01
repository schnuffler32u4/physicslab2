import numpy as np
from functions import M
from scipy.optimize import curve_fit, least_squares
import matplotlib.pyplot as plt

V0 = 10
xdata, ydata = np.loadtxt("fittry.csv", delimiter=",", skiprows=1, unpack=True, encoding='UTF8')
xdata = 1000 * xdata / (2 * np.pi)
ydata = ydata * np.sqrt(2)


def func(x, a, b, c, d, e):
    return 1e10 * np.sqrt((M * x ** 2 * V0 * (d * x - b * x ** 3)) ** 2 / (
            (e - c * x ** 2 + a * x ** 4 - M * x ** 4) ** 2 + (d * x - b * x ** 3) ** 2) ** 2 + (
                                  M * x ** 2 * V0 * (e - c * x ** 2 + a * x ** 4 - M * x ** 4)) ** 2 / (
                                  (e - c * x ** 2 + a * x ** 4 - M * x ** 4) ** 2 - (d * x - b * x ** 3) ** 2) ** 2)


popt, pcov = curve_fit(func, xdata, ydata, p0=[0.8, 1, 5.48, 3.69047095e+13, 1e20], maxfev=500000)
print(M)
print(popt)
# print(pcov)
plt.plot(xdata, func(xdata, *popt), label='fit')
plt.plot(xdata, ydata)
plt.legend()
plt.show()
