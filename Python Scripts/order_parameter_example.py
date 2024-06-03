import matplotlib.pyplot as plt
import numpy as np

from guassian_order import calc_order_param_weighted, random_interval
np.random.seed(23452987)
def get_angles(base_ang, angle_interval, neigh_size, base_len, length_interval):
    angle = random_interval(base_ang - angle_interval / 2, base_ang + angle_interval / 2, neigh_size ** 2)
    lengths = random_interval(base_len - length_interval / 2, base_len + length_interval / 2, neigh_size ** 2)
    return angle, lengths

def get_angles2(angle, lengths):
    return np.cos(np.radians(angle))*lengths, np.sin(np.radians(angle))*lengths

def plot_dist(neigh_size=3, base_ang =45, angle_interval=50, base_len=0.4, length_interval=0.5, fixed_centre=None):
    centre_idx = (neigh_size**2-1)//2
    colors = ['k' for _ in range(neigh_size ** 2)]
    colors[centre_idx] = 'y'

    grid = np.meshgrid(np.linspace(0, 1,  neigh_size), np.linspace(0, 1, neigh_size))
    grid2 = grid[:]

    fig, ax = plt.subplots(dpi=600, figsize=(4.8, 4.8))

    grid1_angles = get_angles(base_ang, angle_interval, neigh_size, base_len, length_interval)
    grid2_angles = get_angles(base_ang, angle_interval, neigh_size, base_len, length_interval)

    if fixed_centre:
        grid1_angles[0][centre_idx] = fixed_centre[0][0]
        grid1_angles[1][centre_idx] = fixed_centre[0][1]
        grid2_angles[0][centre_idx] = fixed_centre[-1][0]
        grid2_angles[1][centre_idx] = fixed_centre[-1][1]

    ax.quiver(*grid, *get_angles2(*grid1_angles), pivot="tail", color=colors, scale=9, headlength=0, headwidth=0, headaxislength=0)
    ax.quiver(*grid2, *get_angles2(*grid2_angles), pivot="tail",color=colors, scale=9, headlength=0, headwidth=0, headaxislength=0)
    ax.scatter(*grid, marker='o', c=colors, s=5)
    ax.scatter(*grid, marker='o', s=1500, facecolors='none', edgecolors='k')

    ax.set_xticks(np.linspace(0, 1,  neigh_size)+1/(neigh_size-1)/2)
    ax.set_yticks(np.linspace(0, 1,  neigh_size)+1/(neigh_size-1)/2)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.tick_params(which="major", bottom=False, left=False)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlim(-0.15, 1.125)
    ax.set_ylim(-0.15, 1.125)
    ax.grid(True)


    centre = []
    all_angles = []
    all_intensities = []

    for _ in range(100):
        angle, lengths = get_angles(base_ang, angle_interval, neigh_size, base_len, length_interval)
        angle = list(angle)
        lengths = list(lengths)


        centre.append([angle.pop(centre_idx), lengths.pop(centre_idx), 0])
        all_angles.append(angle)
        all_intensities.append(lengths)

    ax.set_title(f'S = {calc_order_param_weighted(np.array(all_angles), np.array(all_intensities), fixed_centre if fixed_centre else centre):.1f}')
    fig.tight_layout()
    fig.show()

#plot_dist(angle_interval=0, length_interval=0)
#plot_dist(angle_interval=90, length_interval=0)
plot_dist(base_ang=90, angle_interval=180, length_interval=0.4, fixed_centre=[[100, 0.4, 0], [80,0.4,0]]) # 0
plot_dist(base_ang=45, angle_interval=0, length_interval=0, fixed_centre=[[-40, 0.4, 0], [-50,0.4,0]]) # -1
plot_dist(base_ang=45, angle_interval=0, length_interval=0, fixed_centre=[[40, 0.4, 0], [50,0.4,0]]) # 1
#plot_dist(angle_interval=180, length_interval=0.4)