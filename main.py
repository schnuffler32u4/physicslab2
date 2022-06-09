import numpy as np
from scipy.optimize import curve_fit, least_squares
import matplotlib.pyplot as plt
import functions
from matplotlib.widgets import Slider

plt.style.use('extensys-ms')  # If extensys plots is not installed comment this line

# Importing the data


bcf_d12, bcf_d12_err, bcf1, bcf1_err, bcf2, bcf2_err = np.loadtxt('BCFresults1-2.csv', delimiter=',',
                                                                  skiprows=1, unpack=True)
bcf_d13, bcf_d13_err, bcf13, bcf13_err = np.loadtxt('BCFresults1-3.csv', delimiter=',', skiprows=1,
                                                    unpack=True)
bcf_d3, bcf_d3_err, bcf3, bcf3_err = np.loadtxt('BCFresults3.csv', delimiter=',', skiprows=1, unpack=True)

bcv_d12, bcv_d12_err, bcv1, bcv1_err, bcv2, bcv2_err = np.loadtxt('BCVresults1-2.csv', delimiter=',',
                                                                  skiprows=1, unpack=True)
bcv_d13, bcv_d13_err, bcv13, bcv13_err = np.loadtxt('BCVresults1-3.csv', delimiter=',', skiprows=1,
                                                    unpack=True)
bcv_d3, bcv_d3_err, bcv3, bcv3_err = np.loadtxt('BCVresults3.csv', delimiter=',', skiprows=1, unpack=True)

nrf_d12, nrf_d12_err, nrf1, nrf1_err, nrf2, nrf2_err = np.loadtxt('NRFresults1-2.csv', delimiter=',',
                                                                  skiprows=1, unpack=True)
nrf_d13, nrf_d13_err, nrf13, nrf13_err = np.loadtxt('NRFresults1-3.csv', delimiter=',', skiprows=1,
                                                    unpack=True)
nrf_d3, nrf_d3_err, nrf3, nrf3_err = np.loadtxt('NRFresults3.csv', delimiter=',', skiprows=1, unpack=True)

nrv_d12, nrv_d12_err, nrv1, nrv1_err, nrv2, nrv2_err = np.loadtxt('NRVresults1-2.csv', delimiter=',',
                                                                  skiprows=1, unpack=True)
nrv_d13, nrv_d13_err, nrv13, nrv13_err = np.loadtxt('NRVresults1-3.csv', delimiter=',', skiprows=1,
                                                    unpack=True)
nrv_d3, nrv_d3_err, nrv3, nrv3_err = np.loadtxt('NRVresults3.csv', delimiter=',', skiprows=1, unpack=True)

xdata, xerr, ydata, yerr = np.loadtxt("FM50results.csv", delimiter=",", skiprows=1, unpack=True, encoding='UTF8')
xdata = 1000 * xdata * 2 * np.pi
ydata = ydata * np.sqrt(2)
newx, newxerror, newy, newyerror = np.loadtxt("FM350results.csv", delimiter=",", skiprows=1, unpack=True,
                                              encoding='UTF8')
newx = 1000 * newx * 2 * np.pi
newy = newy * np.sqrt(2)

xfr, _, yfr, _ = np.loadtxt("NRFresults1-3.csv", delimiter=",", skiprows=1, unpack=True, encoding='UTF8')
xfr = xfr / 10
yfr = 1000 * 2 * np.pi * yfr

xdis, ydis = np.loadtxt("vdis.csv", delimiter=",", skiprows=1, unpack=True, encoding='UTF8')
xdis = xdis / 10
ydis = ydis * np.sqrt(2)

# Plotting the data vs the predicted values for a distance of 5 cm
xqr = np.linspace(xdata[0], xdata[-1], 10000)
plt.plot(xqr, functions.voltvar(xqr, 0.1, 0.1, 57, 7, 1e-10), 'r-', label='Prediction')
plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, c='c')
plt.scatter(xdata, ydata, label='Measurements', c='c', s=15, marker="^")
# plt.plot(xdata, functions.voltvar(xdata, 0.0409, 0.086, 57, 7, 1.106e-10), label='Fit but with extracted RLC values')
plt.xlabel('ω [rad/s]')
plt.ylabel('Peak voltage [V]')
plt.legend(loc="upper right", borderaxespad=0.5)
# plt.title("Theoretical model compared with data of frequency versus voltage for a distance of 5 cm")
plt.savefig("Theoretical model compared with data of frequency versus voltage for a distance of 5 cm.jpg", dpi=500)
# plt.show()

# Plotting the data vs the predicted values for a distance of 35 cm
newxqr = np.linspace(newx[0], newx[-1], 10000)
plt.plot(newxqr, functions.voltdis(newxqr, 0.1, 0.1, 57, 7, 1e-10, functions.mutual(0.35)), 'r-', label='Prediction')
plt.errorbar(newx, newy, xerr=newxerror, yerr=newyerror, c='c')
plt.scatter(newx, newy, label='Measurements', c='c', s=15, marker="^")
# plt.plot(xdata, functions.voltvar(xdata, 0.0409, 0.086, 57, 7, 1.106e-10), label='Fit but with extracted RLC values')
plt.xlabel('ω [rad/s]')
plt.ylabel('Peak voltage [V]')
plt.legend(loc="upper right", borderaxespad=0.5)
# plt.title("Theoretical model compared with data of frequency versus voltage for a distance of 35 cm")
plt.savefig("Theoretical model compared with data of frequency versus voltage for a distance of 35 cm.jpg", dpi=500)
# plt.show()

# Fitting the data and plotting with the correct resistance for 5 cm
norespopt, norespcov = curve_fit(functions.voltnores, xdata, ydata, p0=[0.1, 0.1, 1e-10], maxfev=500000)
norespopt = abs(norespopt)  # making sure the values are positive
print(
    f"with L1 = {norespopt[0]}±{functions.roundup(np.sqrt(np.diag(norespcov))[0])} H, L2 = {norespopt[1]}±{functions.roundup(np.sqrt(np.diag(norespcov))[1])} H, C1 = {norespopt[2]}±{functions.roundup(np.sqrt(np.diag(norespcov))[2])} F")
plt.plot(xqr, functions.voltnores(xqr, *norespopt), 'r-', label='Fit with L1 = 0.022±0.004 H, L2 = 0.0913±0.0003 '
                                                                'H, C = 1.4±0.2e-10 F')
plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, c='c')
plt.scatter(xdata, ydata, label='Measurements', c='c', s=15, marker="^")
plt.xlabel('ω [rad/s]')
plt.ylabel('Peak voltage [V]')
plt.legend(loc="upper right", borderaxespad=0.5)
# plt.title("Fitting the values with the correct resistances for a distance of 5 cm")
plt.savefig("Fitting the values with the correct resistances for a distance of 5 cm.jpg", dpi=500)
# plt.show()

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
# plt.title("Fitting the values with the correct resistances for a distance of 35 cm")
plt.savefig("Fitting the values with the correct resistances for a distance of 35 cm.jpg", dpi=500)
# plt.show()

# Fitting the data for all values for 5 cm
print(functions.M)
popt, pcov = curve_fit(functions.voltvar, xdata, ydata, p0=[2.3, 0.1, 57, 7, 1e-10], maxfev=500000)
popt = abs(popt)  # making sure the values are positive
print(
    f"Full fit at 5 cm: L1 = {popt[0]}±{functions.roundup(np.sqrt(np.diag(pcov))[0])} H, L2 = {popt[1]}±{functions.roundup(np.sqrt(np.diag(pcov))[1])} H, R1 = {popt[2]}±{functions.roundup(np.sqrt(np.diag(pcov))[2])} Ω, R2 = {popt[3]}±{functions.roundup(np.sqrt(np.diag(pcov))[3])} Ω, C1 = {popt[4]}±{functions.roundup(np.sqrt(np.diag(pcov))[4])} F")
plt.plot(xqr, functions.voltvar(xqr, *popt), 'r-',
         label='Fit with L1 = 0.0403±0.0002 H, L2 = 0.08715±0.00002 H, R1 = 1067±3 Ω, R2 = 1034±5 Ω, '
               'C1 = 1.128±0.005e-10 F')
plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, c='c')
plt.scatter(xdata, ydata, label='Measurements', c='c', s=15, marker="^")
plt.xlabel('ω [rad/s]')
plt.ylabel('Peak voltage [V]')
plt.legend(loc="upper right", borderaxespad=0.5)
# plt.title("Fitting the values for a distance of 5 cm")
plt.savefig("Fitting the values for a distance of 5 cm.jpg", dpi=500)
# plt.show()
a = popt
# Fitting the data for all values for 35 cm
functions.M = functions.mutual(0.35)
print(functions.M)
popt, pcov = curve_fit(functions.voltvar, newx, newy, p0=[0.1, 0.1, 57, 7, 1e-10], maxfev=500000)
popt = abs(popt)  # making sure the values are positive
print(
    f"Full fit at 35 cm: L1 = {popt[0]}±{functions.roundup(np.sqrt(np.diag(pcov))[0])} H, L2 = {popt[1]}±{functions.roundup(np.sqrt(np.diag(pcov))[1])} H, R1 = {popt[2]}±{functions.roundup(np.sqrt(np.diag(pcov))[2])} Ω, R2 = {popt[3]}±{functions.roundup(np.sqrt(np.diag(pcov))[3])} Ω, C1 = {popt[4]}±{functions.roundup(np.sqrt(np.diag(pcov))[4])} F")
functions.M = functions.mutual(0.35)
plt.plot(newxqr, functions.voltvar(newxqr, *popt), 'r-', label='Fit with ')
plt.errorbar(newx, newy, xerr=newxerror, yerr=newyerror, c='c')
plt.scatter(newx, newy, label='Measurements', c='c', s=15, marker="^")
plt.xlabel('ω [rad/s]')
plt.ylabel('Peak voltage [V]')
plt.legend(loc="upper right", borderaxespad=0.5)
# plt.title("Fitting the values for a distance of 35 cm")
plt.savefig("Fitting the values for a distance of 35 cm.jpg", dpi=500)
# plt.show()

p, = plt.plot(xdata, functions.volt(xdata, *a), 'r-')
plt.ylim(0, 100)
axslider = plt.axes([0.1, 0, 0.8, 0.05])
slider = Slider(ax=axslider, valmin=5, valmax=35, valinit=5, label='Distance (cm)')


def value_update(val):
    M = functions.mutual(slider.val / 100)
    p.set_ydata(functions.voltdis(xdata, *functions.letters(*a), M))


slider.on_changed(value_update)
# plt.show()

"""popt1, pcov1 = curve_fit(functions.v1dis, xdis, ydis, p0=[5, 0.1, 60, 7, 1e-10], bounds=(0, 70), maxfev=5000)
errs = np.sqrt(np.diag(abs(pcov1)))
print(
    f"Full fit at over length: L1 = {popt1[0]}±{functions.roundup(errs[0])} H, L2 = {popt1[1]}±{functions.roundup(errs[1])} H, R1 = {popt1[2]}±{functions.roundup(errs[2])} Ω, R2 = {popt1[3]}±{functions.roundup(errs[3])} Ω, C1 = {popt1[4]}±{functions.roundup(errs[4])} F")

plt.plot(xdis, ydis, c='c', label='data')
q, = plt.plot(xdis, functions.v1dis(xdis, *popt1), 'r-', label='Fit over distance with ')
plt.plot(xdis, functions.v1dis(xdis, 4, 0.1, 57, 7, 1e-10), 'g-', label='Interesting parameters pizda ma-tii de sci-py ')


# bxslider = plt.axes([0.2, 0.1, 0.65, 0.05])
# cxslider = plt.axes([0.2, 0, 0.65, 0.05])
# Lslider = Slider(ax=bxslider, valmin=0, valmax=5, valinit=0.1, label='L1 [H]')
# Rslider = Slider(ax=cxslider, valmin=0, valmax=60, valinit=57, label='R1 [Ω]')


def newupdate(val):
    q.set_ydata(v1dis(xdis, *functions.letters(Lslider.val, 0.1, 57, Rslider.val, 1e-10)))


# Lslider.on_changed(newupdate)
# Rslider.on_changed(newupdate)
plt.xlabel('Distance [cm]')
plt.ylabel('First peak voltage [V]')
plt.legend()
# plt.show()"""

popt1, pcov1 = curve_fit(functions.f1dis, xfr, yfr, p0=[0.04, 0.087, 1067, 1034, 1.128e-10], maxfev=50000)
errs = np.sqrt(np.diag(abs(pcov1)))
errs = [1, 1, 1, 1, 1]
print(
    f"Full fit frequency vs length: L1 = {popt1[0]}±{functions.roundup(errs[0])} H, L2 = {popt1[1]}±{functions.roundup(errs[1])} H, R1 = {popt1[2]}±{functions.roundup(errs[2])} Ω, R2 = {popt1[3]}±{functions.roundup(errs[3])} Ω, C1 = {popt1[4]}±{functions.roundup(errs[4])} F")
plt.plot(xfr, yfr, c='c', label='data')
plt.plot(xfr, functions.f1dis(xfr, *popt1), 'r-', label='Fit over distance with ')
plt.plot(xfr, functions.f1dis(xfr, 0.1, 0.1, 57, 7, 1e-10), 'g-',
         label='Interesting parameters pizda ma-tii de sci-py ')
plt.xlabel('Distance [cm]')
plt.ylabel('First peak angular frequency [rad/s]')
plt.legend()
# plt.show()

dist = np.linspace(50, 350, 100)

functions.double_plot(nrv_d13, nrf_d13, nrv13, nrf13, nrv_d13_err, nrf_d13_err, nrv13_err, nrf13_err, dist,
                      functions.f1dis(dist / 10, 0.01, 0.01, 57, 7, 1e-10) / (2000 * np.pi), dist,
                      functions.v1dis(dist / 10, 0.01, 0.01, 57, 7, 1e-10) / np.sqrt(2), "first, merged", "no resistor")

functions.double_plot(bcv_d13, bcf_d13, bcv13, bcf13, bcv_d13_err, bcf_d13_err, bcv13_err, bcf13_err, dist,
                      functions.f1dis(dist / 10, 0.01, 0.01, 57, 1000, 1e-10) / (2000 * np.pi), dist,
                      functions.v1dis(dist / 10, 0.01, 0.01, 57, 1000, 1e-10) / np.sqrt(2), "first, merged", "resistor")

functions.double_plot(nrv_d12, nrf_d12, nrv2, nrf2, nrv_d12_err, nrf_d12_err, nrv2_err, nrf2_err, dist,
                      functions.f2dis(dist / 10, 0.01, 0.01, 57, 7, 1e-10) / (2000 * np.pi), dist,
                      functions.v2dis(dist / 10, 0.01, 0.01, 57, 7, 1e-10) / np.sqrt(2), "second", "no resistor")

functions.double_plot(bcv_d12, bcf_d12, bcv2, bcf2, bcv_d12_err, bcf_d12_err, bcv2_err, bcf2_err, dist,
                      functions.f1dis(dist / 10, 0.01, 0.01, 57, 1007, 1e-10) / (2000 * np.pi), dist,
                      functions.v1dis(dist / 10, 0.01, 0.01, 57, 1007, 1e-10) / np.sqrt(2), "second", "resistor")
