import h5py
from matplotlib.mlab import bivariate_normal
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors
import pandas as pd

from aux_11 import *
from plotter import *
from parser import define_variables
import helper as hlp

from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

root = '../outputs/run_01/'

compare_nobw = False
if compare_nobw:
    name_nobw = 'gr200'
    root_nobw = '/tigress/PERSEUS/hakobyan/TRISTAN/dec_2017/' + folder + '/' + name_nobw + '/'
    root_nobw += 'output/'

output_dir = root + 'pics/'
# name = 'test_error'
# root = '/tigress/PERSEUS/hakobyan/TRISTAN/dec_2017/' + name + '/'
# name = 'gr70_gc50_s400_ppc10_bw20_cl2000_t-4'
# output_dir = '/tigress/hakobyan/pics/december/' + name + '/'

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

    if compare_nobw:
        photons_nobw = getPhotons(root_nobw, step)
        plasma_nobw = getPlasma(root_nobw, step)

    photons = getPhotons(root, step)
    plasma = getPlasma(root, step)
    pairs = plasma[plasma.ind < 0]
    plasma = plasma[plasma.ind > 0]

    yglob_mid = (plasma.y.max() + plasma.y.min()) * 0.5
    xglob_mid = (plasma.x.max() + plasma.x.min()) * 0.5

    x = (np.arange(len(dens[0])) - max(np.arange(len(dens[0]))) * 0.5) * code_downsampling / skin_depth
    y = (np.arange(len(dens)) - max(np.arange(len(dens))) * 0.5) * code_downsampling / skin_depth
    x, y = np.meshgrid(x, y)

    fig = plt.figure(figsize=(26, 13))
    global_fontsize = 15

    ax00 = plt.subplot2grid((3,6),(0,0),colspan=3)
    ax10 = plt.subplot2grid((3,6),(1,0),colspan=3)
    ax20 = plt.subplot2grid((3,6),(2,0),colspan=3)
    ax01 = plt.subplot2grid((3,6),(0,3),colspan=3)
    ax11 = plt.subplot2grid((3,6),(1,3),colspan=3)
    ax21 = plt.subplot2grid((3,6),(2,3),colspan=3)
    # ax21 = plt.subplot2grid((3,6),(2,3),colspan=3)

    xmin = x.min()
    xmax = x.max()
    ymin = -100 / (skin_depth / 10.)
    ymax = 100 / (skin_depth / 10.)

    ax00.streamplot(x, y, by, bx, color = (1,1,1,0.5), density = 0.6, linewidth = 1.2)

    cbar00, ax00 = plot_dens(ax00, x, y,
                     dens / ppc0, vmin = 1, vmax = 5e3, label = r'plasma $[n_0]$',
                     xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                     cmap = 'plasma',
                     scaling = 'log', setover = 'red', extend = 'both', ret_cbar = True)

#     cnt = ax00.contour(x, y, dens / ppc0, levels = [bw_dens_lim / ppc0], colors = '#a3ff82', linewidths = 1.0)
#     cbar00.add_lines(cnt, erase=False)
#     cnt.collections[0].set_label('p-p limit')
#     if cool_dens_lim < 5e3:
#         cnt = ax00.contour(x, y, dens / ppc0, levels = [cool_dens_lim / ppc0], colors = '#6df2ff', linewidths = 1.0)
#         cbar00.add_lines(cnt, erase=False)
#         (cbar00.ax.get_children()[1]).set_linewidth(2)
#         (cbar00.ax.get_children()[2]).set_linewidth(2)
#         cnt.collections[0].set_label('cooling limit')
#     else:
#         (cbar00.ax.get_children()[1]).set_linewidth(2)
#     legend = ax00.legend(loc='upper right',fontsize=global_fontsize, facecolor='white')
#     legend.get_frame().set_facecolor('#ffffff')
#     legend.get_frame().set_alpha(0.9)
#     for legobj in legend.legendHandles:
#         legobj.set_linewidth(2.0)

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
    legend = ax01.legend(loc='upper right',fontsize=global_fontsize, facecolor='white')
    legend.get_frame().set_facecolor('#ffffff')
    legend.get_frame().set_alpha(0.9)
    (cbar01.ax.get_children()[1]).set_linewidths(2)
    for legobj in legend.legendHandles:
        legobj.set_linewidth(2.0)

#     ax20 = plot_temperature(ax20, plasma,
#                     xmin, xmax, ymin, ymax,
#                     max_g = 50, skin_depth = skin_depth, dwn = 8)

#     exmax = max(ex.max(), -ex.min())**0.3
#     ax20 = plot_dens(ax20, x, y,
#                      np.sign(ex) * np.abs(ex)**0.3, vmin = -exmax, vmax = exmax, label = r'$\pm|E_x|^{0.3}$',
#                      xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
#                      cmap = 'bwr',
#                      scaling = 'lin', extend = 'neither')

    ax11 = plot_dens(ax11, x, y,
                     dnpair, vmin = 1e-2, vmax = 10, label = r'pairs formed at a given timestep',
                     xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                     cmap = 'plasma',
                     scaling = 'log', setover = 'red', extend = 'both')

#     ax11 = plot_stat(ax11, root, step,
#                     epsph_min, epsph_max)

    if photons is not None:
        original = len(photons)
        if original > 0:
            if len(photons) > 1e7:
                photons = photons.sample(1e7)
            ax21 = plot_spectrum(ax21, photons.e, stride = stride * original / len(photons),
                                 label='photons', color = 'black',
                                 weights = photons.ch, max_e = 1e4)
#         if compare_nobw:
#             original = len(photons_nobw)
#             if original > 0:
#                 if len(photons_nobw) > 1e7:
#                     photons_nobw = photons_nobw.sample(1e7)
#                 ax21 = plot_spectrum(ax21, photons_nobw.e, stride = 100 * stride * original / len(photons_nobw),
#                                      label='photons (no pp)', color = 'black', ls = '--',
#                                      weights = photons_nobw.ch, max_e = 1e4)
    # [(plasma.x > xglob_mid - 50 * skin_depth) & (plasma.x < xglob_mid + 50 * skin_depth)]
    if compare_nobw:
        ax21 = plot_spectrum(ax21, plasma_nobw.g, stride=stride,
                             label='plasma (no pp)', color = (0,0,1,1), ls = '--', max_e = 1e4)
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
