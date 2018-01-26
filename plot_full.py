import h5py
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
import pandas as pd

from aux_11 import *
from plotter import *
from parser import define_variables

root = '../../outputs/test_sd3/'

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

start = 0
end = max_number

print "There are overall approx " + str(max_number) + " files."

for step in range(start, min(max_number, end), 5):
    print step
    dens = getField(root, step, 'dens', getSizes(root, step))
    densph = getField(root, step, 'densph', getSizes(root, step))
    denbw = getField(root, step, 'denbw', getSizes(root, step))
    dnpair = getField(root, step, 'dnpair', getSizes(root, step))
    bx = getField(root, step, 'bx', getSizes(root, step))
    by = getField(root, step, 'by', getSizes(root, step))
    bz = getField(root, step, 'bz', getSizes(root, step))
    ex = getField(root, step, 'ex', getSizes(root, step))
    ey = getField(root, step, 'ey', getSizes(root, step))
    ez = getField(root, step, 'ez', getSizes(root, step))
    ecrossb = np.sqrt((ey * bz - ez * by)**2 + (ex * bz - bz * ex)**2 + (ex * by - ey * bx)**2)
    bsquared = bx**2 + by**2 + bz**2
    ecrossb /= bsquared
    multiplicity = divideArray(denbw, dens - denbw)

    photons = getPhotons(root, step)
    plasma = getPlasma(root, step)
    pairs = plasma[plasma.ind < 0]
    plasma = plasma[plasma.ind > 0]

    yglob_mid = (plasma.y.max() + plasma.y.min()) * 0.5
    xglob_mid = (plasma.x.max() + plasma.x.min()) * 0.5

    x = (np.arange(len(dens[0])) - max(np.arange(len(dens[0]))) * 0.5) * code_downsampling / skin_depth
    y = (np.arange(len(dens)) - max(np.arange(len(dens))) * 0.5) * code_downsampling / skin_depth
    x, y = np.meshgrid(x, y)

    fig = plt.figure(figsize=(15, 15))
    global_fontsize = 15

    ax00 = plt.subplot2grid((3,6),(0,0),colspan=3)
    ax10 = plt.subplot2grid((3,6),(1,0),colspan=3)
    ax20 = plt.subplot2grid((3,6),(2,0),colspan=3)
    ax01 = plt.subplot2grid((3,6),(0,3),colspan=3)
    ax11 = plt.subplot2grid((3,6),(1,3),colspan=3)
    ax21 = plt.subplot2grid((3,6),(2,3),colspan=3)

    xmin = x.min()
    xmax = x.max()
    ymin = y.min()
    ymax = y.max()
    # ymin = -100 / (skin_depth / 10.)
    # ymax = 100 / (skin_depth / 10.)

    ax00.streamplot(x, y, by, bx, color = (1,1,1,0.5), density = 0.6, linewidth = 1.2)

    cbar00, ax00 = plot_dens(ax00, x, y,
                     dens / ppc0, vmin = 1, vmax = 5e3, label = r'plasma $[n_0]$',
                     xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                     cmap = 'plasma',
                     scaling = 'log', setover = 'red', extend = 'both', ret_cbar = True)

    ax10 = plot_dens(ax10, x, y,
                     densph, vmin = 10, vmax = 1e5, label = 'photons',
                     xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                     cmap = 'afmhot',
                     scaling = 'log', extend = 'both')

    cbar01, ax01 = plot_dens(ax01, x, y,
                     multiplicity, vmin = 1e-1, vmax = 100, label = 'pair-plasma multiplicity',
                     xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                     cmap = 'inferno',
                     scaling = 'log', setover = 'white', extend = 'both', ret_cbar = True)

    ax20 = plot_dens(ax20, x, y,
                     np.sqrt(bsquared / np.mean(bsquared[1])), vmin = 1e-2, vmax = 30, label = r'$B/B_{\rm up}$',
                     xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                     cmap = 'inferno',
                     scaling = 'log', setover = 'red', extend = 'both')

    cnt = ax01.contour(x, y, multiplicity, levels = [1.], colors = '#f4b942', linewidths = 1.0)
    cbar01.add_lines(cnt, erase=False)
    cnt.collections[0].set_label('multiplcity = 1')
    legend = ax01.legend(loc='upper right',fontsize=global_fontsize)
    legend.get_frame().set_facecolor('#ffffff')
    legend.get_frame().set_alpha(0.9)

    ax11 = plot_dens(ax11, x, y,
                     dnpair, vmin = 1e-2, vmax = 10, label = r'pairs formed at a given timestep',
                     xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                     cmap = 'plasma',
                     scaling = 'log', setover = 'red', extend = 'both')

    if photons is not None:
        original = len(photons)
        if original > 0:
            if len(photons) > 1e7:
                photons = photons.sample(1e7)
            ax21 = plot_spectrum(ax21, photons.e, stride = stride * original / len(photons),
                                 label='photons', color = 'black',
                                 weights = photons.ch, max_e = 1e4)

    if len(pairs) > 0:
        ax21 = plot_spectrum(ax21, pairs.g, stride=stride,
                             label='pairs', color = (1,0,0,1), max_e = 1e4)
        ax21 = plot_spectrum(ax21, np.concatenate((plasma.g, pairs.g)),
                             stride=stride,
                             label='plasma', color = 'blue', max_e = 1e4)
    else:
        ax21 = plot_spectrum(ax21, plasma.g,
                             stride=stride,
                             label='plasma', color = 'blue', max_e = 1e4)

    plt.tight_layout()
    plt.savefig(output_dir + "all_" + str(step+1).zfill(3) + ".png", dpi=150)
