import os
import pickle

import numpy as np
import skimage as ski
from matplotlib import pyplot as plt
from guassian_order import find_first


def get_angles2(angle, lengths):
    return np.cos(np.radians(angle+90))*lengths, np.sin(np.radians(angle+90))*lengths

pickle_path = r"./pickles/DP1_21cm/window300_outmost.pickle"
image_path = r"D:\Data\Microscopy\RCM\2023\20230926 - DP1\21cm\outmost3\outmost3_MMStack_Pos0.ome.tif"
csv_path = r"./input/DP1_21cm/window300_outmost.csv"

# 0 1   2       3       4           5           6       7           8:
# x	y	width	height	Max Index	Mask Median	Angle	Relevance   data:
RFT = np.loadtxt(csv_path, delimiter=',', skiprows=1)

image = ski.io.imread(image_path)

if os.path.isfile(pickle_path):
    with open(pickle_path, 'rb') as pickle_file:
        results = pickle.load(pickle_file)
        print(f"Loaded pickle with {len(results)} rows")

fig, ax = plt.subplots()
ax.imshow(image, cmap='gray')

base_nan = (np.nan, np.nan, np.nan)
match_angle = 35
buff = 0
st = 0 - buff
en = 180 + buff

np.random.seed(23452987)

coords = [RFT[i, :2] for i in range(RFT.shape[0])]
window = RFT[0, 3]
RFT_sum = np.sum(RFT[:, 8:], axis=0)
p_per_a = RFT_sum.shape[0] / 180
pts = int((en - st) * p_per_a)
angles = np.linspace(st, en, pts)
# parameters [angle, intensity, variance]x n, std error, residuals, min_val


xs, ys = [np.unique(l) for l in zip(*coords)]
img_intensity_map = [[[] for _ in range(len(xs))] for _ in range(len(ys))]
peak_angle_map = [[[] for _ in range(len(xs))] for _ in range(len(ys))]

X = []
Y = []
U = []
V = []
for i in range(RFT.shape[0]):
    x = find_first(xs, RFT[i, 0])
    y = find_first(ys, RFT[i, 1])


    # select: sig pars, bg pars, std error, residuals, min_val, image median intensity, main peaks match
    res = list(results[i])
    angle_length = np.sum([lst[1] for lst in res[0] + res[1] if lst])
    angle = res[0][np.array(res[0])[:, 1].argmax()][0] if res[0] else np.nan
    if angle_length > 10:
        (ax, ay) = get_angles2(angle, angle_length)
        X.append(RFT[i, 0] + window // 2)
        Y.append(RFT[i, 1] + window // 2)
        U.append(ax)
        V.append(ay)


plt.quiver(X, Y, U, V, color="y", pivot="middle", headlength=0, headaxislength=0, headwidth=0)

plt.show()