import matplotlib.pyplot as plt

x = ["OJ", "OJ-split", "RFT", "OJ", "OJ-split", "RFT", "AFT", "AFT"]
y = [-6.1, -6.5, -11.4, 8.8, 9.1, 6.7, -10.3, 9]
error = [7.9, 9.3, 9.5, 5.5, 6.9, 7.5, 8, 8.7]

plt.errorbar(x, y, yerr=error, fmt='o', capsize=5)
plt.ylabel("angle ($^\circ$)")
plt.show()