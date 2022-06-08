import numpy as np
from scipy.optimize import curve_fit, least_squares
import matplotlib.pyplot as plt
import functions
from matplotlib.widgets import Slider

plt.style.reload_library()
plt.style.use('extensys-ms')

xdata, xerr, ydata, yerr = np.loadtxt("FM50results.csv", delimiter=",", skiprows=1, unpack=True, encoding='UTF8')
xdata = 1000 * xdata * 2 * np.pi
ydata = ydata * np.sqrt(2)
newx, newxerror, newy, newyerror = np.loadtxt("FM350results.csv", delimiter=",", skiprows=1, unpack=True,
                                              encoding='UTF8')
newx = 1000 * newx * 2 * np.pi
newy = newy * np.sqrt(2)

# Plotting the data vs the measured values for a distance of 5 cm
plt.plot(xdata, functions.voltvar(xdata, 0.1, 0.1, 57, 7, 1e-10), 'r-', label='Prediction')
plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, c='c')
plt.scatter(xdata, ydata, label='Measurements', c='c', s=15, marker="^")
# plt.plot(xdata, functions.voltvar(xdata, 0.0409, 0.086, 57, 7, 1.106e-10), label='Fit but with extracted RLC values')
plt.xlabel('ω [rad/s]')
plt.ylabel('Peak voltage [V]')
plt.legend(loc="upper right", borderaxespad=0.5)
plt.title("Theoretical model compared with data of frequency versus voltage for a distance of 5 cm")
plt.show()

# Plotting the data vs the measured values for a distance of 35 cm
plt.plot(newx, functions.voltdis(newx, 0.1, 0.1, 57, 7, 1e-10, functions.mutual(35)), 'r-', label='Prediction')
plt.errorbar(newx, newy, xerr=newxerror, yerr=newyerror, c='c')
plt.scatter(newx, newy, label='Measurements', c='c', s=15, marker="^")
# plt.plot(xdata, functions.voltvar(xdata, 0.0409, 0.086, 57, 7, 1.106e-10), label='Fit but with extracted RLC values')
plt.xlabel('ω [rad/s]')
plt.ylabel('Peak voltage [V]')
plt.legend(loc="upper right", borderaxespad=0.5)
plt.title("Theoretical model compared with data of frequency versus voltage for a distance of 35 cm")
plt.show()

# Fitting the data and plotting with the correct resistance for 5 cm
norespopt, norespcov = curve_fit(functions.voltnores, xdata, ydata, p0=[0.1, 0.1, 1e-10], maxfev=500000)
norespopt = abs(norespopt)  # making sure the values are positive
print(
    f"with L1 = {norespopt[0]}±{functions.roundup(np.sqrt(np.diag(norespcov))[0])} H, L2 = {norespopt[1]}±{functions.roundup(np.sqrt(np.diag(norespcov))[1])} H, C1 = {norespopt[2]}±{functions.roundup(np.sqrt(np.diag(norespcov))[2])} F")
plt.plot(xdata, functions.voltnores(xdata, *norespopt), 'r-', label='Fit with L1 = 0.022±0.004 H, L2 = 0.0913±0.0003 '
                                                                    'H, C = 1.4±0.2e-10 F')
plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, c='c')
plt.scatter(xdata, ydata, label='Measurements', c='c', s=15, marker="^")
plt.xlabel('ω [rad/s]')
plt.ylabel('Peak voltage [V]')
plt.legend(loc="upper right", borderaxespad=0.5)
plt.show()

# Fitting the data and plotting with the correct resistance for 35 cm
norespopt, norespcov = curve_fit(functions.voltnores35, newx, newy, p0=[0.1, 0.1, 1e-10], maxfev=500000)
norespopt = abs(norespopt)  # making sure the values are positive
print(
    f"with L1 = {norespopt[0]}±{functions.roundup(np.sqrt(np.diag(norespcov))[0])} H, L2 = {norespopt[1]}±{functions.roundup(np.sqrt(np.diag(norespcov))[1])} H, C1 = {norespopt[2]}±{functions.roundup(np.sqrt(np.diag(norespcov))[2])} F")
plt.plot(newx, functions.voltnores35(newx, *norespopt), 'r-',
         label='Fit with L1 = 0±10000 H, L2 = 0±400 H, C1 = 0±8e-6')
plt.errorbar(newx, newy, xerr=newxerror, yerr=newyerror, c='c')
plt.scatter(newx, newy, label='Measurements', c='c', s=15, marker="^")
plt.xlabel('ω [rad/s]')
plt.ylabel('Peak voltage [V]')
plt.legend(loc="upper right", borderaxespad=0.5)
plt.show()

# Fitting the data for all values for 5 cm
print(functions.M)
popt, pcov = curve_fit(functions.voltvar, xdata, ydata, p0=[2.3, 0.1, 57, 7, 1e-10], maxfev=500000)
popt = abs(popt)  # making sure the values are positive
print(
    f"Full fit at 5 cm: L1 = {popt[0]}±{functions.roundup(np.sqrt(np.diag(pcov))[0])} H, L2 = {popt[1]}±{functions.roundup(np.sqrt(np.diag(pcov))[1])} H, R1 = {popt[2]}±{functions.roundup(np.sqrt(np.diag(pcov))[2])} Ω, R2 = {popt[3]}±{functions.roundup(np.sqrt(np.diag(pcov))[3])} Ω, C1 = {popt[4]}±{functions.roundup(np.sqrt(np.diag(pcov))[4])} F")
plt.plot(xdata, functions.voltvar(xdata, *popt), 'r-',
         label='Fit with L1 = 0.0403±0.0002 H, L2 = 0.08715±0.00002 H, R1 = 1067±3 Ω, R2 = 1034±5 Ω, C1 = 1.128±0.005e-10 F')
plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, c='c')
plt.scatter(xdata, ydata, label='Measurements', c='c', s=15, marker="^")
plt.xlabel('ω [rad/s]')
plt.ylabel('Peak voltage [V]')
plt.legend(loc="upper right", borderaxespad=0.5)
plt.show()

# Fitting the data for 35 cm
functions.M = functions.mutual(35) * 100
print(functions.M)
popt, pcov = curve_fit(functions.voltvar, newx, newy, p0=[0.1, 0.1, 57, 7, 1e-10], maxfev=500000)
popt = abs(popt)  # making sure the values are positive
print(
    f"Full fit at 35 cm: L1 = {popt[0]}±{functions.roundup(np.sqrt(np.diag(pcov))[0])} H, L2 = {popt[1]}±{functions.roundup(np.sqrt(np.diag(pcov))[1])} H, R1 = {popt[2]}±{functions.roundup(np.sqrt(np.diag(pcov))[2])} Ω, R2 = {popt[3]}±{functions.roundup(np.sqrt(np.diag(pcov))[3])} Ω, C1 = {popt[4]}±{functions.roundup(np.sqrt(np.diag(pcov))[4])} F")
functions.M = functions.mutual(35)
plt.plot(newx, functions.voltvar(newx, *popt), 'r-', label='Fit with ')
plt.errorbar(newx, newy, xerr=newxerror, yerr=newyerror, c='c')
plt.scatter(newx, newy, label='Measurements', c='c', s=15, marker="^")
plt.xlabel('ω [rad/s]')
plt.ylabel('Peak voltage [V]')
plt.legend(loc="upper right", borderaxespad=0.5)
plt.show()

p, = plt.plot(xdata, functions.volt(xdata, *popt), 'r-')
plt.ylim(0, 100)
axslider = plt.axes([0.1, 0, 0.8, 0.05])
slider = Slider(ax=axslider, valmin=5, valmax=35, valinit=5, label='Distance (cm)')


def value_update(val):
    M = functions.mutual(slider.val / 100)
    p.set_ydata(functions.voltdis(xdata, *functions.letters(*popt), M))


slider.on_changed(value_update)
plt.show()

xdis, ydis = np.loadtxt("vdis.csv", delimiter=",", skiprows=1, unpack=True, encoding='UTF8')
xdis = xdis / 10
ydis = ydis * np.sqrt(2)


def vdis(x, L1, L2, R1, R2, C1):
    xfreq = np.linspace(0, 1000000, 1000)
    y = np.zeros(len(x))
    yfreq = np.zeros(len(x))
    for i in range(len(x)):
        y[i] = max(functions.voltdis(xfreq, *functions.letters(L1, L2, R1, R2, C1), functions.mutual(x[i] / 100)))
        for r in functions.voltdis(xfreq, *functions.letters(L1, L2, R1, R2, C1), functions.mutual(x[i] / 100)):
            if r == y[i]:
                yfreq[i] = xfreq[i]
    return y


def fdis(x, L1, L2, R1, R2, C1):
    xfreq = np.linspace(0, 1000000, 1000)
    y = np.zeros(len(x))
    yfreq = np.zeros(len(x))
    for i in range(len(x)):
        y[i] = max(functions.voltdis(xfreq, *functions.letters(L1, L2, R1, R2, C1), functions.mutual(x[i] / 100)))
        for r in functions.voltdis(xfreq, *functions.letters(L1, L2, R1, R2, C1), functions.mutual(x[i] / 100)):
            if r == y[i]:
                yfreq[i] = xfreq[i]
    return yfreq


popt1, pcov1 = curve_fit(vdis, xdis, ydis, p0=[5, 0.1, 60, 7, 1e-10], bounds=(0,70), maxfev=5000)
errs = np.sqrt(np.diag(abs(pcov1)))
print(
    f"Full fit at over length: L1 = {popt1[0]}±{functions.roundup(errs[0])} H, L2 = {popt1[1]}±{functions.roundup(errs[1])} H, R1 = {popt1[2]}±{functions.roundup(errs[2])} Ω, R2 = {popt1[3]}±{functions.roundup(errs[3])} Ω, C1 = {popt1[4]}±{functions.roundup(errs[4])} F")

plt.plot(xdis, ydis, c='c', label='data')
# plt.plot(xdis, vdis(xdis, *popt1), label='fit')
q, = plt.plot(xdis, vdis(xdis, *popt1), 'r-', label='Fit over distance with ')
plt.plot(xdis, vdis(xdis, 4, 0.1, 57, 7, 1e-10), 'g-', label='Interesting parameters pizda ma-tii de sci-py ')

# bxslider = plt.axes([0.2, 0.1, 0.65, 0.05])
# cxslider = plt.axes([0.2, 0, 0.65, 0.05])
# Lslider = Slider(ax=bxslider, valmin=0, valmax=5, valinit=0.1, label='L1 [H]')
# Rslider = Slider(ax=cxslider, valmin=0, valmax=60, valinit=57, label='R1 [Ω]')


def newupdate(val):
    q.set_ydata(vdis(xdis, *functions.letters(Lslider.val, 0.1, 57, Rslider.val, 1e-10)))


# Lslider.on_changed(newupdate)
# Rslider.on_changed(newupdate)
plt.xlabel('Distance [cm]')
plt.ylabel('First peak voltage [V]')
plt.legend()
plt.show()
