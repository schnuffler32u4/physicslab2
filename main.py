import numpy as np
from scipy.optimize import curve_fit, least_squares
import matplotlib.pyplot as plt
import functions
from matplotlib.widgets import Slider
# plt.style.use('extensys')

V0 = 10.3
xdata, ydata = np.loadtxt("fittry.csv", delimiter=",", skiprows=1, unpack=True, encoding='UTF8')
xdata = 1000 * xdata * 2 * np.pi
ydata = ydata * np.sqrt(2)
newx, newy = np.loadtxt("fitcheck.csv", delimiter=",", skiprows=1, unpack=True, encoding='UTF8')
newx = 1000 * newx * 2 * np.pi
newy = newy * np.sqrt(2)

popt, pcov = curve_fit(functions.volt, xdata, ydata, p0=[0.01, 1.4, 2e9, 1.4e14, 1e20], maxfev=500000)
print(functions.M)

print(popt)
print(functions.circuit_values(*popt))
print(functions.letters(*functions.circuit_values(*popt)))
plt.plot(xdata, functions.volt(xdata, *popt), label='Fit')
plt.plot(xdata, ydata, label='DZata')
plt.plot(xdata, functions.voltvar(xdata, *functions.circuit_values(*popt)), label='Fit but with extracted RLC values')
plt.xlabel('ω [rad/s]')
plt.ylabel('Peak voltage [V]')
plt.legend()
plt.show()

plt.plot(newx, newy)
plt.plot(newx, functions.voltdis(newx, *popt, functions.mutual(0.35)), label='fit')
plt.legend()
plt.show()

p, = plt.plot(xdata, functions.volt(xdata, *popt))
plt.ylim(0, 100)
axslider = plt.axes([0.1, 0, 0.8, 0.05])
slider = Slider(ax=axslider, valmin=5, valmax=35, valinit=5, label='Distance (cm)')


def value_update(val):
    M = functions.mutual(slider.val / 100)
    p.set_ydata(functions.voltdis(xdata, *popt, M))


slider.on_changed(value_update)
plt.show()

xdis, ydis = np.loadtxt("vdis.csv", delimiter=",", skiprows=1, unpack=True, encoding='UTF8')
ydis = ydis * np.sqrt(2)


def vdis(x, a, b, c, d, e):
    xfreq = np.linspace(0, 1000000, 10000)
    y = np.zeros(len(x))
    for i in range(len(x)):
        y[i] = max(functions.voltdis(xfreq, a, b, c, d, e, functions.mutual(x[i] / 100)))
    return y


popt1, pcov1 = curve_fit(vdis, xdis, ydis, p0=[*functions.letters(0.1, 0.1, 60, 7, 1e-10)], maxfev=500000)
print(functions.circuit_values(*popt1))
plt.plot(xdis, ydis)
plt.plot(xdis, vdis(xdis, *popt1), label='fit')

plt.legend()
plt.show()
