import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from guassian_order import flatten

header = ["Window", "width" ,"angle", "angle", "mean ang", "std ang" ,"mean int", "std int", "mean std", "std std"]

data = np.loadtxt('characteristic_size.csv', skiprows=1)

FOV = 512 # px
windows = np.unique(data[:,0])
wavenumbers = np.unique(data[:,1])
wavelengths = wavenumbers*FOV*2 # px
angles = np.unique(data[:,3])
mean_std = np.zeros((windows.shape[0],wavelengths.shape[0]))
std_std = np.zeros_like(mean_std)
fit_values = []
x_fit = []
for x, window in enumerate(windows):
    for y, wavenumber in enumerate(wavenumbers):
        core_mean_std = data[np.logical_and(data[:,0]==window,data[:,1]==wavenumber)][:, 8]
        core_std_std = data[np.logical_and(data[:,0]==window,data[:,1]==wavenumber)][:, 9]
        mean_std[x,y] = np.average(core_mean_std)
        std_std[x,y] = np.sqrt(np.var(core_mean_std)**2+np.average(core_std_std)**2)

for x, window in enumerate(windows):
    x_vals = window/wavelengths
    relevant_x = np.logical_and(x_vals > 0.5, x_vals < 5) # nyquist
    plt.errorbar(x_vals, mean_std[x,:], yerr=std_std[x, :], label=window, fmt="o", markersize=5)
    x_fit += list(x_vals[relevant_x])
    fit_values += list(mean_std[x, relevant_x])

ax = plt.gca()
ax.set_xscale('log')



params = [-15.3, 26.8]
x_indices = np.logspace(-1, 0.75, 100)

#fit_eq = lambda x, *params: params[1]*np.exp(x*params[0]) # a = -0.65, b = 50.4
fit_eq = lambda x, *params: params[0]*np.log(x)+params[1]
params, covariance = curve_fit(fit_eq, x_fit, fit_values, p0=params)
points = fit_eq(x_indices, *params)
print(*params)
plt.plot(x_indices, points, label="fit", zorder=10)
ax.vlines(0.5, 0, 1, transform=ax.get_xaxis_transform(), colors='r', label="Nyquist")
plt.legend()
plt.ylabel("σ(°)")
plt.xlabel("waves per window j ($px^{-1}$)")
plt.show()