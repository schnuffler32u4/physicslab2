import numpy as np
from scipy.optimize import curve_fit, least_squares
import matplotlib.pyplot as plt
import functions
from matplotlib.widgets import Slider

V0 = 10.3
xdata, ydata = np.loadtxt("fittry.csv", delimiter=",", skiprows=1, unpack=True, encoding='UTF8')
xdata = 1000 * xdata * 2 * np.pi
ydata = ydata * np.sqrt(2)
newx, newy = np.loadtxt("fitcheck.csv", delimiter=",", skiprows=1, unpack=True, encoding='UTF8')
newx = 1000 * newx * 2 * np.pi
newy = newy * np.sqrt(2)
M = functions.mutual(0.05)


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


def sanity(L1, L2, R1, R2, C1):
    a = L1 * L2
    b = L1 * R2 + R1 * L2
    c = L1 / 1e-10 + R1 * R2 + L2 / C1
    d = R1 / 1e-10 + R2 / C1
    e = 1 / (C1 * 1e-10)
    return [a, b, c, d, e]


def pula(x, a, b, c, d, e, M):
    s = a * x ** 4 - M ** 2 * x ** 4 - c * x ** 2 + e
    t = d * x - b * x ** 3
    p = M * x ** 2 * V0
    A = p * t / (s ** 2 + t ** 2)
    B = p * s / (s ** 2 + t ** 2)
    return 1e10 * np.sqrt(A ** 2 + B ** 2)


popt, pcov = curve_fit(func, xdata, ydata, p0=[0.01, 1.4, 2e9, 1.4e14, 1e20], maxfev=500000)
print(M)

print(popt)
print(functions.circuit_values(*popt))
print(sanity(*functions.circuit_values(*popt)))
plt.plot(xdata, func(xdata, *popt), label='fit')
plt.plot(xdata, ydata)
plt.plot(xdata, funcvar(xdata, *functions.circuit_values(*popt)), label='yummy')
plt.legend()
plt.show()

# plt.plot(newx, newy)
# plt.plot(newx, func(newx, 0.01, 1.4, 2e9, 1.4e11, 1e20), label='fit')
# plt.legend()
# plt.show()

p, = plt.plot(xdata, func(xdata, *popt))
axslider = plt.axes([0.1, 0, 0.8, 0.05])
slider = Slider(ax=axslider, valmin=5, valmax=35, valinit=5, label='Distance (cm)')


def value_update(val):
    M = functions.mutual(slider.val / 100)
    p.set_ydata(pula(xdata, *popt, M))


slider.on_changed(value_update)
plt.show()
