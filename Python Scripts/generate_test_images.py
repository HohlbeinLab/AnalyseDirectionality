import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib as mpl

from guassian_order import calc_order_param_weighted, flatten

image_folder = "./test_images_2"

# scroll through angles and widths

data_type = 16 # bpp of generates images
resolution = 512 # width/height of the image
tiling = 4 # increase if artefacts appear

min_width = 0.01 # fraction of resolution
max_width = 1 # fraction of resolution
delta_width = 0.05

angle_delta = 5 # deg

def generate_sequence():
    for width in np.arange(min_width, max_width, delta_width): # from large to small
        x_values = np.array(256 * (0.5 + 0.5 * np.cos((np.pi * (x_indices - resolution)) / (0.5 * width * resolution))), dtype=np.dtype("i"))
        img = np.tile(x_values, (2*resolution, 1))
        image = Image.fromarray(img)

        for angle in np.arange(0, 45+angle_delta, angle_delta): # from vertical to horizontal to vertical
            rot_image = image.rotate(angle)
            crop_image = rot_image.crop((resolution//2, resolution//2, resolution+resolution//2, resolution+resolution//2))
            crop_image.save(f"{image_folder}/width{width:.2f}_angle{angle}.tif")

def generate_order_image():
    widths = [[0.11, 0.11, 0.11], [0.11, 0.11, 0.11], [0.11, 0.11, 0.11]]
    #angles = [[-1, -1, -1], [-1, -1, -1], [-1, -1, -1]]
    #angles = [[130, 130, 130], [130, 130, 130], [130, 130, 130]]
    #angles = [[130, 130, 130], [130, 40, 130], [130, 130, 130]]
    angles = [[40, 40, 40], [130, 130, 130], [130, 130, 130]]
    og_angles = [l[:] for l in angles]
    intensities =  [[1, 1, 1], [1, 1, 1], [1, 1, 1]]
    #intensities = [[1, 1, 1], [1, 1, 1], [1, 1, 1]]
    size_divisor = len(angles)
    data = []
    r = resolution // size_divisor
    x_indices = np.arange(0, tiling * resolution, 1, dtype=np.dtype("i"))  # 0, 1, 2, ..., resolution
    WOPs = []

    for y in range(size_divisor):
        for x in range(len(angles[0])):
            if angles[y][x] == -1:
                angles[y][x] = np.random.random()*180
            x_values = np.array(
                intensities[y][x]* 256 * (0.5 + 0.5 * np.cos((np.pi * (x_indices - resolution)) / (0.5 * widths[y][x] * resolution))),
                dtype=np.dtype("i"))
            img = np.tile(x_values, (tiling * resolution, 1))
            image = Image.fromarray(img)
            rot_image = image.rotate(angles[y][x])
            crop_image = rot_image.crop(((x+0.5+tiling/2)*r, (y+0.5+tiling/2)*r, (x+1.5+tiling/2)*r, (y+1.5+tiling/2)*r))
            data.append(crop_image)

    total_image = np.concatenate([
        np.concatenate(data[:size_divisor], axis=1),
        np.concatenate(data[size_divisor:size_divisor*2], axis=1),
        np.concatenate(data[-3:], axis=1)
    ], axis=0)
    plt.imshow(total_image, cmap=mpl.colormaps.get_cmap('gray'))

    centre = [
        [og_angles[size_divisor // 2].pop(size_divisor // 2), intensities[size_divisor // 2].pop(size_divisor // 2),
         widths[size_divisor // 2].pop(size_divisor // 2)]]

    for _ in range(10000):
        angles = [[b if b != -1 else np.random.random()*180 for b in a] for a in og_angles]

        if centre[0][0] == -1:
            centre[0][0] = np.random.random() * 180

        WOP = calc_order_param_weighted(flatten(angles), flatten(intensities), centre)
        WOPs.append(WOP)

    ax = plt.gca()
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.tick_params(which="major", bottom=False, left=False)
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    print(np.average(WOPs))
    print(np.var(WOPs))
    plt.title(f'{np.average(WOPs):.2f}±{np.var(WOPs):.2f}')
    plt.savefig(rf"H:\PhD\AnalyseDirectionality\Python Scripts\Figures\WOP\varied_{np.average(WOPs):.2f}.svg")
    plt.show()

generate_order_image()