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

import sys
sys.path.append('/u/hhakoby1/vis/tqdm/')
from tqdm import tqdm

root = '/u/hhakoby1/outputs/production_runs/s3e3_gc150/'
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

for step in tqdm(range(start, min(max_number, end) + 1, 5)):
    dens = getField(root, step, 'dens', getSizes(root, step), ymin = 30, ymax = -30)
    densph = getField(root, step, 'densph', getSizes(root, step), ymin = 30, ymax = -30)
    denbw = getField(root, step, 'denbw', getSizes(root, step), ymin = 30, ymax = -30)
    dnpair = getField(root, step, 'dnpair', getSizes(root, step), ymin = 30, ymax = -30)
    bx = getField(root, step, 'bx', getSizes(root, step), ymin = 30, ymax = -30)
    by = getField(root, step, 'by', getSizes(root, step), ymin = 30, ymax = -30)
    bz = getField(root, step, 'bz', getSizes(root, step), ymin = 30, ymax = -30)
    bsquared = bx**2 + by**2 + bz**2
    multiplicity = divideArray(denbw, dens - denbw)

    x = (np.arange(len(dens[0])) - max(np.arange(len(dens[0]))) * 0.5) * code_downsampling / skin_depth
    y = (np.arange(len(dens)) - max(np.arange(len(dens))) * 0.5) * code_downsampling / skin_depth
    x, y = np.meshgrid(x, y)

    fig = plt.figure(figsize=(20, 16))
    global_fontsize = 15

    ax00 = plt.subplot2grid((3,4),(0,0),colspan=2)
    ax10 = plt.subplot2grid((3,4),(1,0),colspan=2)
    ax20 = plt.subplot2grid((3,4),(2,0),colspan=2)
    ax01 = plt.subplot2grid((3,4),(0,2),colspan=2)
    ax11 = plt.subplot2grid((3,4),(1,2),colspan=2)
    ax210 = plt.subplot2grid((3,4),(2,2))
    ax211 = plt.subplot2grid((3,4),(2,3))

    xmin = x.min()
    xmax = x.max()
    ymin = -200 / (skin_depth / 10.)
    ymax = 200 / (skin_depth / 10.)

    ax00.streamplot(x, y, by, bx, color = (1,1,1,0.5), density = 0.6, linewidth = 1.2)

    cbar00, ax00 = plot_dens(ax00, x, y,
                     dens / ppc0, vmin = 1, vmax = 5e3, label = r'plasma $[n_0]$',
                     xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                     cmap = 'plasma',
                     scaling = 'log', setover = 'red', extend = 'both', ret_cbar = True)

    ax10 = plot_dens(ax10, x, y,
                     densph, vmin = 10, vmax = 1e7, label = 'photons',
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

    ax11 = plot_dens(ax11, x, y,
                     dnpair, vmin = 1e-2, vmax = 100, label = r'pairs formed at a given timestep',
                     xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                     cmap = 'plasma',
                     scaling = 'log', setover = 'red', extend = 'both')

    data = h5py.File(root + 'spec.tot.{}'.format(str(step).zfill(3)), 'r')
    bins = data['bn'].value
    parts = data['npart'].value
    pairs = data['npair'].value
    phots = data['nphot'].value

    param = h5py.File(root + 'param.{}'.format(str(step).zfill(3)), 'r')
    nprocs = param['sizey'].value[0] * param['sizex'].value[0]

    max_e = 1e5
    min_n = 1e0
    max_n = 1e10
    ax211 = plot_spectrum_new(ax211, bins, phots, nprocs,
                              label = 'photons', color = 'black',
                              max_e = max_e, min_n = min_n, max_n = max_n)

    min_e = 1
    max_e = 1e5
    min_n = 1e0
    max_n = 1e8
    ax210 = plot_spectrum_new(ax210, bins, parts, nprocs,
                              label = 'all plasma', color = 'blue',
                              min_e = min_e, max_e = max_e, min_n = min_n, max_n = max_n)
    ax210 = plot_spectrum_new(ax210, bins, pairs, nprocs,
                              label = 'pairs', color = 'red',
                              min_e = min_e, max_e = max_e, min_n = min_n, max_n = max_n)

    plt.tight_layout()
    plt.savefig(output_dir + "all_" + str(step).zfill(3) + ".png", dpi=150)
