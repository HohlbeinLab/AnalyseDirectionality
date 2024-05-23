import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit


def gauss(x, *params):
   return params[1] * np.exp(-((x - params[0]) / params[2]) ** 2)


def func(x, *params):
    y = np.zeros_like(x)

    y += params[0] * np.exp(-params[1] * x)
    if len(params) > 2:
        y += params[2] * np.exp(-params[3] * x)

    return y

RCM = np.loadtxt("SPC_140_coronal_size_2", delimiter='\t', skiprows=0)
RCM[:, 1] = RCM[:, 1] - 0.3

plt.plot(RCM[:, 0], RCM[:, 1], c='b', label = 'data')

guess = [
    25, -1 / 25,

]

guess_gauss = [
    7, 1, 0.7,
    10, 2, 10,

]


popt, pcov = curve_fit(func, RCM[:, 0], RCM[:, 1], p0=guess)



pred = func(RCM[:, 0], *popt)
func2 = lambda x, *params: pred + sum([(params[i+1] * np.exp(-((x - params[i]) / params[i+2]) ** 2)) for i in range(0, len(params), 3)])
popt2, pcov2 = curve_fit(func2, RCM[:, 0], RCM[:, 1], p0=guess_gauss)

plt.plot(RCM[:, 0], func2(RCM[:, 0], *popt2), 'g-', label='sum fit')
plt.plot(RCM[:, 0], gauss(RCM[:, 0], *popt2[0:3]), 'r-', label='gauss fit')
plt.plot(RCM[:, 0], gauss(RCM[:, 0], *popt2[3:6]), 'r-')
plt.plot(RCM[:, 0], func(RCM[:, 0], *popt), 'r-')

print(popt)
print(popt2)


plt.ylim([-1, 5])
plt.grid(True)
plt.legend()
plt.ylabel("intensity (a.u.)")
plt.xlabel("wavenumber (1/mm)")
plt.show()

print(f"size: {1/popt2[0]:.2f}mm +- {abs(1000/(popt2[2]+popt2[0])-(1000/popt2[0])):.2f}um ")
print(f"size: {1/popt2[3]:.2f}mm +- {abs(1000/(popt2[2]+popt2[0])-(1000/popt2[0])):.2f}um ")