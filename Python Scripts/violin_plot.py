import numpy as np
import matplotlib.pyplot as plt
import os
from guassian_order import weighted_nanmean
import matplotlib.patches as mpatches
import seaborn as sns

labels=[]


loadpaths = ["./results/DP1/DP1_3cm", "./results/DP1/DP1_21cm"]
label = ["3 cm", "21 cm"]
sides = ["high", "high"]
offset = [0, 0.25]
colors = ["blue", "red"]
positions2 = [[6, 5, 4, 3, 2, 1], [6, 2, 3, 4, 5, 1]]
#positions2 = [[7, 8, 9, 10, 11, 12], [7, 8, 9, 10, 11, 12]]
for i, loadpath in enumerate(loadpaths):
    data = []
    positions = np.array(positions2[i])
    weighted_means = []
    for filename in os.listdir(loadpath):
        file_path = fr"{loadpath}\{filename}\raw\{filename}_order_all angles_neighbourhood_9_raw.csv"
        weights_path = fr"{loadpath}\{filename}\raw\{filename}_peak_intensity_raw.csv"
        A = np.loadtxt(file_path, delimiter=',').flatten()
        weights = np.loadtxt(weights_path, delimiter=',').flatten()
        filter = np.logical_and(~np.isnan(A), weights > np.nanmean(weights)*0)
        filtered_A = A[filter]
        if not filtered_A.shape[0]:
            filtered_A = [0]
        data.append(filtered_A)
        weighted_means.append(weighted_nanmean(filtered_A, weights[filter]))

    parts = plt.violinplot(data, positions=1*(positions+offset[i]), showextrema=False, side=sides[i], showmeans=False)
    labels.append((mpatches.Patch(color=colors[i]), label[i]))
    #sns.swarmplot(data=data, size=2)

    for pc in parts['bodies']:
        pc.set_facecolor(colors[i])
        pc.set_edgecolor('black')
        pc.set_alpha(0.6)

    quartile3 = [np.percentile(d, [75]) for d in data]
    median = [np.percentile(d, [50]) for d in data]
    plt.scatter(positions+offset[i]+0.1, weighted_means, marker='o', color=colors[i], s=30, zorder=3)
    #plt.scatter(positions+offset[i]+0.1, quartile3, marker='o', color=colors[i], s=30, zorder=3)

plt.legend(*zip(*labels))
plt.ylabel("WOP (-)")
plt.xlabel("SP")
#plt.ylim((0, 0.7))
plt.show()
plt.savefig("./figures/A2_violin_plot.svg")