import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def func(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i + 1]
        wid = params[i + 2]
        y = y + amp * np.exp(-((x - ctr) / wid) ** 2)
    return y


def plot_fit(ax, popt, arr, ind, comment, col='r-'):
    for i in range(0, len(popt), 3):
        print(f"{comment}: {popt[i]:.1f}+-{popt[i+2]:.1f} degrees ({popt[i+1]:.1f})")
        line, = ax.plot(arr[:ind], func(arr[:ind], *popt[i:i+3]), col)
    return line

def MRI(arr):

    fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=False, sharex=True)
    lim = 30
    ax1.set_xlim(-lim, lim)
    ax2.set_xlim(-lim, lim)
    ax2.set_xlabel("angle ($^\circ$)")
    ax1.set_ylabel("counts")
    ax2.set_ylabel("counts")


    counts, bins = np.histogram(arr[:, 2], bins=180)

    guess = [
        -10, 40, 10,
        0, 20, 30,
        10, 40, 10
    ]

    popt, pcov = curve_fit(func, bins[:-1], counts, p0=guess)

    plot_fit(ax1, popt, bins, -1, "sum")


    ax1.plot(bins[:-1], counts, c='b', label='data')

    ax1.plot(bins[:-1], func(bins[:-1], *popt), 'g-', label="sum fit")
    y_cutoff = 128 / 2
    arr1 = arr[arr[:, 1] < y_cutoff]
    arr2 = arr[arr[:, 1] > y_cutoff]

    counts, bins = np.histogram(arr1[:, 2], bins=180)
    ax2.plot(bins[:-1], counts, c='b')
    popt, pcov = curve_fit(func, bins[:-1], counts, p0=[-10, 40, 10])
    line = plot_fit(ax2, popt, bins, -1, "top")
    line.set_label('gauss fit')

    counts, bins = np.histogram(arr2[:, 2], bins=180)
    ax2.plot(bins[:-1], counts, c='b')
    popt, pcov = curve_fit(func, bins[:-1], counts, p0=[10, 40, 10])
    plot_fit(ax2, popt, bins, -1, "bottom")



    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    fig.legend(lines, labels)

    fig.show()



def MRI_own(arr, offset, correction, angle_offset, guess):
    fig, ax = plt.subplots(1, 1, constrained_layout=True)
    start = offset
    arr[:, 1] = arr[:, 1] - correction
    arr[:, 0] = (arr[:, 0] + start)*(360/400)-angle_offset
    ax.plot(arr[:, 0], arr[:, 1], c='b', label='data')

    ax.set_xlabel("angle ($^\circ$)")
    ax.set_ylabel("intensity (a.u)")

    # x center, y height, error


    popt, pcov = curve_fit(func, arr[:, 0], arr[:, 1], p0=guess)


    for i in range(0, len(popt), 3):
        print(f"{popt[i]:.1f}+-{popt[i+2]:.1f} degrees ({popt[i+1]:.1f})")
        line, = ax.plot(arr[:, 0], func(arr[:, 0], *popt[i:i+3]), 'r-')

    line.set_label('gauss fit')
    ax.plot(arr[:, 0], func(arr[:, 0], *popt), 'g-', label='sum fit')

    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    fig.legend(lines, labels)

    fig.show()


def RCM_plot(arr):
    fig, ax = plt.subplots(1, 1)

    counts, bins = np.histogram(arr[:, 5], bins=180)

    guess = [
        0, 80, 25,
    ]

    popt, pcov = curve_fit(func, bins[:-1], counts, p0=guess)

    for i in range(0, len(popt), 3):
        print(f"{popt[i]:.1f}+-{popt[i + 2]:.1f} degrees ({popt[i + 1]:.1f})")
        ax.plot(bins[:-1], func(bins[:-1], *popt[i:i + 3]), 'r-')

    ax.plot(bins[:-1], counts, c='b')
    ax.plot(bins[:-1], func(bins[:-1], *popt), 'g-', label='sum fit')

    fig.show()


def MRI_own2(arr):
    fig, ax = plt.subplots(1, 1, constrained_layout=True)

    arr[:, 1] = arr[:, 1] - 0.9
    arr[:, 0] = (arr[:, 0] ) * (360 / 400) - 90
    ax.plot(arr[:, 0], arr[:, 1], c='b', label='data')

    ax.set_xlabel("angle ($^\circ$)")
    ax.set_ylabel("intensity (a.u)")

    # x center, y height, error
    # Sam image
    guess = [
        0, 0.2, 25,
        15, 0.4, 10,
    ]

    popt, pcov = curve_fit(func, arr[:, 0], arr[:, 1], p0=guess)

    for i in range(0, len(popt), 3):
        print(f"{popt[i]:.1f}+-{popt[i + 2]:.1f} degrees ({popt[i + 1]:.1f})")
        line, = ax.plot(arr[:, 0], func(arr[:, 0], *popt[i:i + 3]), 'r-')

    line.set_label('gauss fit')
    ax.plot(arr[:, 0], func(arr[:, 0], *popt), 'g-', label='sum fit')

    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    fig.legend(lines, labels)

    fig.show()


# X, Y, Orientation, Co. x En.
OJ = np.loadtxt("SPC_140C_coronal_OrientationJ.txt", delimiter=',', skiprows=0)

# 0, 1, 2,     3,  4,  5,           6          7
# X, Y, Slice, Dy, Dy, Orientation, Coherency, Energy
OJ_test = np.loadtxt("test.txt", delimiter='\t', skiprows=0)
Mart = np.loadtxt("SPC_140_coronal_Martijn.txt", delimiter='\t', skiprows=0)
RCM = np.loadtxt("test_RCM", delimiter='\t', skiprows=0)
RCM_2 = np.loadtxt("RCM_3_2", delimiter='\t', skiprows=0)
RFT_airyscan = np.loadtxt("Airyscan_RFT", delimiter='\t', skiprows=0)

RCM = np.c_[RCM, RCM[:, 6]*RCM[:, 7]]
OJ_test = np.c_[OJ_test, OJ_test[:, 6]*OJ_test[:, 7]]
OJ = OJ_test[:, (0, 1, 5, 8)]


#MRI(OJ)

# Sam image


#MRI_own(Mart, 43, 0.35, 90, [-9,  0.3, 5, 0, 0.55, 25,9.34, 0.6,  5,89.62483385,  0.22754284, 25.18526383])
#MRI_own2(RCM_2)
RCM_plot(RCM)

#MRI_own(RFT_airyscan, 130, 0.13, 180, [-25, 0.06, 10, 40, 0.14, 10])

