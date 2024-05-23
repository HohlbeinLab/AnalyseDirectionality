import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

def gauss(x, *params):
   return params[1] * np.exp(-((x - params[0]) / params[2]) ** 2)

def func(x, *params):
    y = np.zeros_like(x)

    y += params[0] * np.exp(-params[1] * x)
    y += params[2] * np.exp(-params[3] * x)

    return y



RCM = np.loadtxt("RCM_3_2_size", delimiter='\t', skiprows=0)
RCM[:, 1] = RCM[:, 1] - 0.9


plt.plot(RCM[:, 0], RCM[:, 1], c='b', label='data')

guess = [
    25, -1 / 25,
    12, -1 / 40,
]

guess_gauss = [10, 5, 0.7]

popt, pcov = curve_fit(func, RCM[:, 0], RCM[:, 1], p0=guess)
pred = func(RCM[:, 0], *popt)
func2 = lambda x, *params: pred + params[1] * np.exp(-((x - params[0]) / params[2]) ** 2)
popt2, pcov2 = curve_fit(func2, RCM[:, 0], RCM[:, 1], p0=guess_gauss)

plt.plot(RCM[:, 0], func2(RCM[:, 0], *popt2), 'g-', label='sum fit')
plt.plot(RCM[:, 0], func(RCM[:, 0], *popt), 'r-', label='fit')
plt.plot(RCM[:, 0], gauss(RCM[:, 0], *popt2), 'r-')
print(popt)
print(popt2)

plt.yscale("log")
plt.ylim( (pow(10,-2),pow(10,1.5)) )
plt.grid(True)

plt.ylabel("intensity (a.u.)")
plt.xlabel("wavenumber (1/um)")
plt.legend()
plt.show()

print(f"size: {1/popt2[0]:.2f}um +- {abs(1000/(popt2[2]+popt2[0])-(1000/popt2[0])):.2f}nm ")