import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



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
    relevant_x = x_vals > 0 # nyquist
    plt.errorbar(x_vals, mean_std[x,:], yerr=std_std[x, :], label=f"{window/FOV:.2f}", fmt="o", markersize=5)
    x_fit += list(x_vals[relevant_x])
    fit_values += list(mean_std[x, relevant_x])


ax = plt.gca()
ax.set_xscale('log')



#params = [-15.3, 26.8]
#params = [-0.65, 50.4]
params = [19.32, 1.95, 1.218, 1.044]
x_indices = np.logspace(-1.3, 1.5, 100)

fit_eq = lambda x, *params: params[3]+params[0]/(1+np.power(x/params[2], params[1]))
invert_fit_eq = lambda x, *params: params[2]*np.power((params[0]/(x-params[3]))-1, 1/params[1])



#fit_eq = lambda x, *params: params[1]*np.exp(x*params[0]) # a = -0.65, b = 50.4
#fit_eq = lambda x, *params: params[0]*np.log(x)+params[1]
params, covariance = curve_fit(fit_eq, x_fit, fit_values, p0=params)
points = fit_eq(x_indices, *params)
print(*params)
plt.plot(x_indices, points, label="fit", zorder=10)
#ax.vlines(0.5, 0, 1, transform=ax.get_xaxis_transform(), colors='r', label="Nyquist")
ax.vlines(1, 0, 1, transform=ax.get_xaxis_transform(), colors='b', label="Nyquist Frequency")
plt.legend()
plt.ylabel("σ(°)")
plt.xlabel("waves per window j ($px^{-1}$)")
plt.savefig("./characteristic_size.svg")
plt.show()

plt.figure()
plt.xlabel("original wavelength (px/FOV)")
plt.ylabel("calculated wavelength (px/FOV)")

for x, window in enumerate(windows):
    plt.scatter(wavelengths/FOV, window/invert_fit_eq(mean_std[x, :], *params)/FOV, label=f"{window}, {window/FOV:.2f}")

plt.ylim([0, 1])
plt.xlim([0, 1])
plt.plot(wavelengths/FOV, wavelengths/FOV, label="unity", zorder=10)
plt.legend()
plt.savefig("./inverse.svg")
plt.show()

input_angles = data[:,2]
detected_angles = data[:,3]
detected_angles[detected_angles > 178] -= 180
plt.figure()
plt.scatter(input_angles, detected_angles)
plt.xlabel("input angles(°)")
plt.ylabel("calculated angles (°)")
plt.savefig("./input_angles.svg")
plt.show()
