import numpy as np
from scipy.optimize import curve_fit, least_squares
import matplotlib.pyplot as plt
import functions
from matplotlib.widgets import Slider

plt.style.reload_library()
plt.style.use('extensys-ms')

V0 = 10.3
xdata, xerr, ydata, yerr = np.loadtxt("FM50results.csv", delimiter=",", skiprows=1, unpack=True, encoding='UTF8')
xdata = 1000 * xdata * 2 * np.pi
ydata = ydata * np.sqrt(2)
newx, newy = np.loadtxt("fitcheck.csv", delimiter=",", skiprows=1, unpack=True, encoding='UTF8')
newx = 1000 * newx * 2 * np.pi
newy = newy * np.sqrt(2)

# popt, pcov = curve_fit(functions.volt, xdata, ydata, p0=[0.01, 1.4, 2e9, 1.4e14, 1e20], maxfev=500000)
print(functions.M)
popt, pcov = curve_fit(functions.voltvar, xdata, ydata, p0=[2.3, 0.1, 57, 7, 1e-10], maxfev=500000)
popt = abs(popt)
print(popt)
# print(functions.circuit_values(*popt))
# print(functions.letters(*functions.circuit_values(*popt)))
print(np.sqrt(np.diag(pcov)))
# plt.plot(xdata, functions.volt(xdata, *popt), label='Fit')
plt.plot(xdata, functions.voltvar(xdata, *popt), label='Fit')
plt.plot(xdata, ydata, label='Data')
# plt.plot(xdata, functions.voltvar(xdata, 0.0409, 0.086, 57, 7, 1.106e-10), label='Fit but with extracted RLC values')
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
xdis = xdis / 10
ydis = ydis * np.sqrt(2)


def vdis(x, a, b, c, d, e):
    xfreq = np.linspace(0, 1000000, 1000)
    y = np.zeros(len(x))
    for i in range(len(x)):
        y[i] = max(functions.voltdis(xfreq, a, b, c, d, e, functions.mutual(x[i] / 100)))
    return y


# popt1, pcov1 = curve_fit(vdis, xdis, ydis, p0=[*functions.letters(0.1, 5, 60, 7, 1e-10)], maxfev=5000)
# print(functions.circuit_values(*popt1))
plt.plot(xdis, ydis, label='data')
# plt.plot(xdis, vdis(xdis, *popt1), label='fit')
q, = plt.plot(xdis, vdis(xdis, *functions.letters(2.7, 0.1, 60, 7, 1e-10)), label='model')
bxslider = plt.axes([0.2, 0.1, 0.65, 0.05])
cxslider = plt.axes([0.2, 0, 0.65, 0.05])
Lslider = Slider(ax=bxslider, valmin=0, valmax=5, valinit=0.1, label='L1 [H]')
Rslider = Slider(ax=cxslider, valmin=0, valmax=60, valinit=57, label='R1 [Ω]')


def newupdate(val):
    q.set_ydata(vdis(xdis, *functions.letters(Lslider.val, 0.1, Rslider.val, 7, 1e-10)))


Lslider.on_changed(newupdate)
Rslider.on_changed(newupdate)
plt.xlabel('Distance [cm]')
plt.ylabel('First peak voltage [V]')
plt.legend()
plt.show()
