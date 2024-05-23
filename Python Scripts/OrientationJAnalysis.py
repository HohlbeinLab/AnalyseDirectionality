import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit


def func(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i + 1]
        wid = params[i + 2]
        y = y + amp * np.exp(-((x - ctr) / wid) ** 2)
    return y

OJ = np.loadtxt("Airyscan.txt", delimiter='\t', skiprows=0)
OJ = np.c_[OJ, OJ[:, 6]*OJ[:, 7]]
OJ = OJ[:, (0, 1, 5, 8)]

bins = []

start = np.min(OJ[:, 2])
end = np.max(OJ[:, 2])
step = 1

angles = np.linspace(start, end, int((end-start)/step))

for angle in angles:
    in_bin = OJ[OJ[:, 2] < angle + step]
    in_bin = in_bin[in_bin[:, 2] > angle]
    bins.append(np.sum(in_bin[:, 3]))

plt.plot(angles, bins, c='b', label='data')

guess = [
    50, 2.5, 10,
]

popt, pcov = curve_fit(func, angles, bins, p0=guess)

for i in range(0, len(popt), 3):
    print(f"fit: {popt[i]:.1f}+-{popt[i+2]:.1f} degrees ({popt[i+1]:.1f})")
    plt.plot(angles, func(angles, *popt[i:i+3]), c='r', label='fit')

plt.xlabel("angle ($^\circ$)")
plt.ylabel("intensity (a.u)")
plt.legend()
plt.show()