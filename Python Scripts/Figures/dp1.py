import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("DP1", skiprows=1)

positions = [3, 21]
for pos in positions:
    plt.errorbar(data[data[:,0]==pos,2], data[data[:,0]==pos,3], yerr=data[data[:,0]==pos,4], fmt="o", markersize=5, label=f"{pos} cm")

plt.legend()
plt.ylabel("mean WOP")
plt.xlabel("Distance from core (μm)")
plt.show()
plt.savefig("meanWOP_DP1.svg")