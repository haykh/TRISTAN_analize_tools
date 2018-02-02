import matplotlib
matplotlib.use("agg")
import h5py
from matplotlib.mlab import bivariate_normal
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors

import sys
sys.path.append('/u/hhakoby1/vis/tqdm/')
from tqdm import tqdm

from aux_11 import *
from plotter import *
from color_data import plasma_cmap
plt.register_cmap(name='plasma', cmap=plasma_cmap)
from parser import define_variables

root = '/u/hhakoby1/outputs/test_big/'
output_dir = root + 'pics/'
simulation_variables = define_variables(root)

import os
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

max_number = int(getNumberOfFiles(root))

code_downsampling = simulation_variables['code_downsampling']
skin_depth = simulation_variables['skin_depth']
ppc0 = simulation_variables['ppc0']
speed_of_light = simulation_variables['speed_of_light']
output_period = simulation_variables['output_period']
bw_dens_lim = simulation_variables['bw_dens_lim']
cool_dens_lim = simulation_variables['cool_dens_lim']
epsph_max = simulation_variables['epsph_max']
epsph_min = simulation_variables['epsph_min']
stride = simulation_variables['stride']

start = 1
end = max_number

for step in tqdm(range(start, min(max_number, end))):
    dens = getField(root, step, 'dens', getSizes(root, step))
    densph = getField(root, step, 'densph', getSizes(root, step))
    denbw = getField(root, step, 'denbw', getSizes(root, step))
    multiplicity = divideArray(denbw, dens - denbw)

    x = (np.arange(len(dens[0])) - max(np.arange(len(dens[0]))) * 0.5) * code_downsampling / skin_depth
    y = (np.arange(len(dens)) - max(np.arange(len(dens))) * 0.5) * code_downsampling / skin_depth
    x,y = np.meshgrid(x,y)
    xmin = x.min()
    xmax = x.max()
    ymin = y.min()
    ymax = y.max()

    fig = plt.figure(figsize=(25,15))

    ax1 = plt.subplot2grid((3,3),(0,0),colspan=3)
    ax2 = plt.subplot2grid((3,3),(1,0),colspan=3)
    ax3 = plt.subplot2grid((3,3),(2,0),colspan=3)

    ax1 = plot_dens(ax1, x, y,
                    dens / ppc0, vmin = 1e-1, vmax = 5e3, label = r'plasma $[n_0]$',
                    xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                    cmap = 'plasma',
                    scaling = 'log', setover = 'red', extend = 'both')
    ax2 = plot_dens(ax2, x, y,
                    densph, vmin = 10, vmax = 1e5, label = 'photons',
                    xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                    cmap = 'afmhot',
                    scaling = 'log', extend = 'both')
    ax3 = plot_dens(ax3, x, y,
                    multiplicity, vmin = 1e-1, vmax = 100, label = 'pair-plasma multiplicity',
                    xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                    cmap = 'Vega10',
                    scaling = 'log', setover = 'white', extend = 'both',

    # plt.tight_layout()
    # plt.savefig(output_dir + "all_" + str(step).zfill(3) + ".png", dpi = 100)
