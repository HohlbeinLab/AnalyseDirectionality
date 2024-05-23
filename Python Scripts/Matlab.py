import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import h5py


def func(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i + 1]
        wid = params[i + 2]
        y = y + amp * np.exp(-((x - ctr) / wid) ** 2)
    return y


AFT = np.loadtxt("results.csv", delimiter=',', skiprows=0)
AFT[:, 2] = (AFT[:, 2] * (180/np.pi)) - 90
AFT = np.c_[AFT, AFT[:, 3]*AFT[:, 4]]


bins = []
simple_bins = []

start = np.min(AFT[:, 2])
end = np.max(AFT[:, 2])
step = 1

angles = np.linspace(start, end, int((end-start)/step))

for angle in angles:
    in_bin = AFT[AFT[:, 2] < angle + step]
    in_bin = in_bin[in_bin[:, 2] > angle]
    bins.append(np.sum(in_bin[:, 5]))
    simple_bins.append(in_bin.shape[0])

plt.plot(angles, simple_bins, c='b', label='data')

guess = [
    50, 2.5, 10,

]

popt, pcov = curve_fit(func, angles, simple_bins, p0=guess)

for i in range(0, len(popt), 3):
    print(f"fit: {popt[i]:.1f}+-{popt[i+2]:.1f} degrees ({popt[i+1]:.1f})")
    plt.plot(angles, func(angles, *popt[i:i+3]), c='r', label='fit')

plt.xlabel("angle ($^\circ$)")
plt.ylabel("counts")
plt.legend()
plt.show()