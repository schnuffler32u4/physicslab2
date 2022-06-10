import scipy.integrate as si
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

V0 = 10.3


def mutual(dis):
    """The function takes the distance in meters as an input and outputs the mutual inductance between the coils."""
    R = 0.105

    def f(x, y):
        return 200 ** 2 * 10 ** (-7) * R ** 2 * np.cos(x - y) / (np.sqrt(2 * R ** 2 * (1 - np.cos(x - y)) + dis ** 2))

    return si.dblquad(f, 0, 2 * np.pi, 0, 2 * np.pi)[0]


M = mutual(0.05)


def volt(x, a, b, c, d, e):
    """Returns the voltage for the given a,b,c,d,e for a value of omega (x), for the M previously defined"""
    s = a * x ** 4 - M ** 2 * x ** 4 - c * x ** 2 + e
    t = d * x - b * x ** 3
    p = M * x ** 2 * V0
    A = p * t / (s ** 2 + t ** 2)
    B = p * s / (s ** 2 + t ** 2)
    return 1e10 * np.sqrt(A ** 2 + B ** 2)


def voltvar(x, L1, L2, R1, R2, C1):
    """Returns the value of the voltage for given circuit values and omega (x), assuming C2 = 1e-10"""
    a = L1 * L2
    b = L1 * R2 + R1 * L2
    c = L1 / 1e-10 + R1 * R2 + L2 / C1
    d = R1 / 1e-10 + R2 / C1
    e = 1 / (C1 * 1e-10)
    return volt(x, a, b, c, d, e)


def letters(L1, L2, R1, R2, C1):
    """Returns the a,b,c,d,e for the given values of the circuits"""
    a = L1 * L2
    b = L1 * R2 + R1 * L2
    c = L1 / 1e-10 + R1 * R2 + L2 / C1
    d = R1 / 1e-10 + R2 / C1
    e = 1 / (C1 * 1e-10)
    return [a, b, c, d, e]


def voltdis(x, a, b, c, d, e, M):
    """Returns the voltage for an omega (x) for a,b,c,d,e variables and a given M (which corresponds to a specific
    distance """
    s = a * x ** 4 - M ** 2 * x ** 4 - c * x ** 2 + e
    t = d * x - b * x ** 3
    p = M * x ** 2 * V0
    A = p * t / (s ** 2 + t ** 2)
    B = p * s / (s ** 2 + t ** 2)
    return 1e10 * np.sqrt(A ** 2 + B ** 2)


def voltdisres(x, a, b, c, d, e, M):
    """Returns the voltage for an omega (x) for a,b,c,d,e variables and a given M (which corresponds to a specific
    distance """
    s = a * x ** 4 - M ** 2 * x ** 4 - c * x ** 2 + e
    t = d * x - b * x ** 3
    p = M * x ** 2 * V0
    A = p * t / (s ** 2 + t ** 2)
    B = p * s / (s ** 2 + t ** 2)
    return 1e3 * x * np.sqrt(A ** 2 + B ** 2)


def circuit_values(a, b, c, d, e):
    """The function takes as input the values of a,b,c,d,e and returns the corresponding capacitances, resistances
    and so on. We assume C2 = 1e-10. The output is an array of the form [L1,L2,R1,R2,C1]"""
    C2 = 1e-10
    C1 = 1 / (e * C2)

    def leq(x):
        return a / (x * C2) + C1 * (b - a * C1 * d / x) / (x - a * C1 / (x * C2)) * (
                d - (b - a * C1 * d / x) / (C2 * (x - a * C1 / (x * C2)))) + x / C1 - c

    L2 = fsolve(leq, 1, maxfev=500000)[0]
    L1 = a / L2
    R1 = (b - a * C1 * d / L2) / (L2 - a * C1 / (L2 * C2))
    R2 = C1 * (d - R1 / C2)
    return [L1, L2, R1, R2, C1]


def voltnores(x, L1, L2, C1):
    """Returns the value of the voltage for given circuit values and omega (x), assuming C2 = 1e-10
    R1 = 57 ohms and R2 = 7 ohms"""
    R1 = 57
    R2 = 7
    a = L1 * L2
    b = L1 * R2 + R1 * L2
    c = L1 / 1e-10 + R1 * R2 + L2 / C1
    d = R1 / 1e-10 + R2 / C1
    e = 1 / (C1 * 1e-10)
    return volt(x, a, b, c, d, e)


def voltnores35(x, L1, L2, C1):
    """Returns the value of the voltage for given circuit values and omega (x), assuming C2 = 1e-10
    R1 = 57 ohms and R2 = 7 ohms"""
    R1 = 57
    R2 = 7
    a = L1 * L2
    b = L1 * R2 + R1 * L2
    c = L1 / 1e-10 + R1 * R2 + L2 / C1
    d = R1 / 1e-10 + R2 / C1
    e = 1 / (C1 * 1e-10)
    return voltdis(x, a, b, c, d, e, mutual(0.35))


def roundup(x):
    """Returns the input value rounded up to one significant figure."""
    if int(np.log10(x)) == np.log10(x):
        y = np.log10(x)
    elif np.log10(x) < 0:
        y = int(np.log10(x)) - 1
    else:
        y = int(np.log10(x))

    if int(x * (10 ** (-y))) * 10 ** y != x:
        return int(x * (10 ** (-y)) + 1) * 10 ** y
    else:
        return x


def v1dis(x, L1, L2, R1, R2, C1):
    """Takes as input the distance and circuit values and returns the voltage of the first peak"""
    xfreq = np.linspace(246e3, 1200e3, 10000)
    y = np.zeros(len(x))
    for i in range(len(x)):
        o = voltdis(xfreq, *letters(L1, L2, R1, R2, C1), mutual(x[i] / 100))
        for j in range(1, len(xfreq) - 1):
            if o[j + 1] < o[j] and o[j - 1] < o[j]:
                y[i] = o[j]
                break
    return y


def f1dis(x, L1, L2, R1, R2, C1):
    """Takes as input the distance and circuit values and returns the frequency of the first peak"""
    xfreq = np.linspace(246e3, 1200e3, 10000)
    y = np.zeros(len(x))
    yfreq = np.zeros(len(x))
    for i in range(len(x)):
        o = voltdis(xfreq, *letters(L1, L2, R1, R2, C1), mutual(x[i] / 100))
        for j in range(1, len(xfreq) - 1):
            if o[j + 1] < o[j] and o[j - 1] < o[j]:
                yfreq[i] = xfreq[j]
                break
    return yfreq


def v2dis(x, L1, L2, R1, R2, C1):
    """Takes as input the distance and circuit values and returns the voltage of the sceond peak"""
    xfreq = np.linspace(246e3, 1200e3, 10000)
    y = np.zeros(len(x))
    for i in range(len(x)):
        o = voltdis(xfreq, *letters(L1, L2, R1, R2, C1), mutual(x[i] / 100))
        for j in range(1, len(xfreq) - 1):
            if o[j + 1] < o[j] and o[j - 1] < o[j]:
                y[i] = o[j]
    return y


def f2dis(x, L1, L2, R1, R2, C1):
    """Takes as input the distance and circuit values and returns the frequency of the second peak"""
    xfreq = np.linspace(246e3, 1200e3, 10000)
    y = np.zeros(len(x))
    yfreq = np.zeros(len(x))
    for i in range(len(x)):
        o = voltdis(xfreq, *letters(L1, L2, R1, R2, C1), mutual(x[i] / 100))
        for j in range(1, len(xfreq) - 1):
            if o[j + 1] < o[j] and o[j - 1] < o[j]:
                yfreq[i] = xfreq[j]
    return yfreq


def v1disres(x, L1, L2, R1, R2, C1):
    """Takes as input the distance and circuit values and returns the voltage of the first peak"""
    xfreq = np.linspace(246e3, 1200e3, 10000)
    y = np.zeros(len(x))
    for i in range(len(x)):
        o = voltdisres(xfreq, *letters(L1, L2, R1, R2, C1), mutual(x[i] / 100))
        for j in range(1, len(xfreq) - 1):
            if o[j + 1] < o[j] and o[j - 1] < o[j]:
                y[i] = o[j]
                break
    return y


def f1disres(x, L1, L2, R1, R2, C1):
    """Takes as input the distance and circuit values and returns the frequency of the first peak"""
    xfreq = np.linspace(246e3, 1200e3, 10000)
    y = np.zeros(len(x))
    yfreq = np.zeros(len(x))
    for i in range(len(x)):
        o = voltdisres(xfreq, *letters(L1, L2, R1, R2, C1), mutual(x[i] / 100))
        for j in range(1, len(xfreq) - 1):
            if o[j + 1] < o[j] and o[j - 1] < o[j]:
                yfreq[i] = xfreq[j]
                break
    return yfreq


def v2disres(x, L1, L2, R1, R2, C1):
    """Takes as input the distance and circuit values and returns the voltage of the sceond peak"""
    xfreq = np.linspace(246e3, 1200e3, 10000)
    y = np.zeros(len(x))
    for i in range(len(x)):
        o = voltdisres(xfreq, *letters(L1, L2, R1, R2, C1), mutual(x[i] / 100))
        for j in range(1, len(xfreq) - 1):
            if o[j + 1] < o[j] and o[j - 1] < o[j]:
                y[i] = o[j]
    return y


def f2disres(x, L1, L2, R1, R2, C1):
    """Takes as input the distance and circuit values and returns the frequency of the second peak"""
    xfreq = np.linspace(246e3, 1200e3, 10000)
    y = np.zeros(len(x))
    yfreq = np.zeros(len(x))
    for i in range(len(x)):
        o = voltdisres(xfreq, *letters(L1, L2, R1, R2, C1), mutual(x[i] / 100))
        for j in range(1, len(xfreq) - 1):
            if o[j + 1] < o[j] and o[j - 1] < o[j]:
                yfreq[i] = xfreq[j]
    return yfreq


def double_plot(x1, x2, y1, y2, x1err, x2err, y1err, y2err, xfmodel, yfmodel, xvmodel, yvmodel, title, resistor):
    fig = plt.figure()

    plt.xticks(fontsize=14)
    plt.xlabel("Coil separation [mm]", fontsize=17)
    # plt.title("Coil separation v " + title + " peak voltage and resonance frequency, " + resistor, fontsize=15)

    ax1 = plt.gca()
    plt.scatter(x2, y2, c='r', label="Resonance frequency", s=23, marker="^")
    plt.errorbar(x2, y2, xerr=x2err, yerr=y2err, c='r')
    plt.plot(xfmodel, yfmodel, 'r-', label="Model prediction for resonance frequency")
    plt.ylabel("AC frequency [kHz]", fontsize=17, color='r')
    plt.yticks(fontsize=14, color='r')

    ax2 = ax1.twinx()
    plt.scatter(x1, y1, c='g', label="Peak rms voltage", s=23, marker="^")
    plt.errorbar(x1, y1, xerr=x1err, yerr=y1err, c='g')
    plt.plot(xvmodel, yvmodel, 'g-', label="Model prediction for rms voltage")
    plt.ylabel("Average rms voltage [V]", fontsize=17, color='g')
    plt.yticks(fontsize=14, color='g')

    fig.legend(loc="upper right", bbox_to_anchor=(0.99, 0.75), bbox_transform=ax1.transAxes, fontsize=15)
    plt.tight_layout()
    # plt.savefig("Coil separation v " + title + " peak voltage and resonance frequency, " + resistor + ".jpg", dpi=500)
    plt.show()

