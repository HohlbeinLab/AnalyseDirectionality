import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from guassian_order import angular_average_weighted, angular_weighted_nanstd

synthetic_path = r'C:\Users\gobes001\LocalSoftware\AnalyseDirectionality\Python Scripts\results\SyntheticRerun'
pickle_path = "./angles.pickle"
angles = defaultdict(list)

if os.path.exists(pickle_path):
    with open(pickle_path, 'rb') as pickle_file:
        angles = pickle.load(pickle_file)
else:
    for filename in os.listdir(synthetic_path):
        if 'xlsx' in filename:
            continue

        angle = int(filename.split('_')[-1])
        np_data = []
        print(filename)
        #angle	 mean ang	 std ang	 mean int	 std int	 mean std	 std std
        csv_path = f"{synthetic_path}/{filename}/csv/{filename}_all_stats.csv"
        data = open(csv_path, 'r').read().split("\n\n")[0].split("\n")[1:]
        for line in data:
            np_data.append(np.fromstring(line, sep=', '))

        np_data = np.array(np_data)
        main_angle = np_data[np_data[:, 3].argmax()]

        angles[angle].append(main_angle[0])
    with open(pickle_path, 'wb') as pickle_file:
        pickle.dump(angles, pickle_file)

for key, value in angles.items():
    print(f'{key}: {angular_average_weighted(value)}')


means = [angular_average_weighted(values) if angular_average_weighted(values) + 2.5 < 180 else angular_average_weighted(values) - 180 for values in angles.values()]
errors = [angular_weighted_nanstd(values) for values in angles.values()]
plt.errorbar(angles.keys(), means, errors, marker="o", linestyle="")
plt.xlabel("input angles (°)")
plt.ylabel("output angles (°)")
plt.plot(angles.keys(), angles.keys(), linestyle="--", label="unity", c="k", zorder=10)
plt.legend()
plt.show()
plt.savefig("./figures/angles_unity.svg")