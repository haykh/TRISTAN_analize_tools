# plot cool_IC

import matplotlib
matplotlib.use("agg")

import sys
sys.path.append('/u/hhakoby1/vis/tqdm/')
from tqdm import tqdm

import matplotlib.pyplot as plt
import numpy as np
import copy
import matplotlib as mpl

from aux_11 import *
from plotter import *
from parser import define_variables
import helper as hlp

# from matplotlib import rc
#
# rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
# rc('text', usetex=True)

root = '/u/hhakoby1/outputs/cool_ic/'
output_dir = root + 'pics/'

simulation_variables = define_variables(root)

import os
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

max_number = int(getNumberOfFiles(root))

code_downsampling = simulation_variables['code_downsampling']
skin_depth = 0.25
ppc0 = simulation_variables['ppc0']
speed_of_light = simulation_variables['speed_of_light']
stride = simulation_variables['stride']

start = 82
end = max_number

global_fontsize = 15

print "There are overall approx " + str(max_number) + " files."

# %matplotlib notebook

for step in tqdm(range(start, min(max_number, end), 1)):
    dens = getField(root, step, 'dens', getSizes(root, step), ymin = 30, ymax = -30)
    densph = getField(root, step, 'densph', getSizes(root, step), ymin = 30, ymax = -30)
    denbw = getField(root, step, 'denbw', getSizes(root, step), ymin = 30, ymax = -30)
    dnpair = getField(root, step, 'dnpair', getSizes(root, step), ymin = 30, ymax = -30)
    bx = getField(root, step, 'bx', getSizes(root, step), ymin = 30, ymax = -30)
    by = getField(root, step, 'by', getSizes(root, step), ymin = 30, ymax = -30)
    bz = getField(root, step, 'bz', getSizes(root, step), ymin = 30, ymax = -30)
    bsquared = bx**2 + by**2 + bz**2
    multiplicity = divideArray(denbw, dens - denbw)

    plasma = getPlasma(root, step)

    yglob_mid = (plasma.y.max() + plasma.y.min()) * 0.5
    xglob_mid = (plasma.x.max() + plasma.x.min()) * 0.5

    dens = np.rot90(dens)

    x = (np.arange(len(dens[0])) - max(np.arange(len(dens[0]))) * 0.5) * code_downsampling / skin_depth
    y = (np.arange(len(dens)) - max(np.arange(len(dens))) * 0.5) * code_downsampling / skin_depth
    x, y = np.meshgrid(x, y)

    fig = plt.figure(figsize=(28, 14))
    global_fontsize = 25

    ax1 = plt.subplot2grid((3,4),(0,0),rowspan=3)
    ax2 = plt.subplot2grid((3,4),(0,1),rowspan=3)
    ax3 = plt.subplot2grid((3,4),(0,2),rowspan=3,colspan=2)

    xmin = -6000
    xmax = 6000
    ymin = y.min()
    ymax = y.max()

    bx = np.rot90(bx)
    by = np.rot90(by)

    # ax1.streamplot(x, y, -bx, by, color = (1,1,1,0.5), density = 2, linewidth = 0.5)

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("top", size="2%", pad=0.05)
    my_cmap = copy.copy(mpl.cm.get_cmap('plasma'))
    my_cmap.set_bad(my_cmap(0))
    my_cmap.set_under(my_cmap(0))
    my_cmap.set_over('red')
    strm = ax1.pcolormesh(x, y, dens, cmap=my_cmap, norm=mpl.colors.LogNorm(vmin=0.1, vmax=1e3))
    cbar = plt.colorbar(strm, cax = cax, extend='both', orientation='horizontal')
    cbar.ax.yaxis.set_tick_params(pad=20)
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.tick_params(labelsize=global_fontsize)
    ax1.tick_params(axis='both', labelsize=global_fontsize)
    ax1.set_aspect(1)
    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(ymin, ymax)
    ax1.set_xlabel(r'$x$, [$c/\omega_{pl}$]', fontsize=global_fontsize)
    ax1.set_ylabel(r'$y$, [$c/\omega_{pl}$]', fontsize=global_fontsize)
    props = dict(boxstyle='square', facecolor='white', alpha=0.9, edgecolor='none')
    ax1.text(0.05, 0.97, r'plasma $[n_0]$', transform=ax1.transAxes, fontsize=global_fontsize, verticalalignment='top', bbox=props)
    ax1.set_xticks([-4000,0,4000])

    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("top", size="2%", pad=0.05)
    my_cmap = copy.copy(mpl.cm.get_cmap('viridis'))
    my_cmap.set_bad(my_cmap(0))
    my_cmap.set_under(my_cmap(0))
    my_cmap.set_over('red')
    strm = ax2.pcolormesh(x, y, np.rot90(np.sqrt(bsquared / np.mean(bsquared[1]))), cmap=my_cmap, norm=mpl.colors.LogNorm(vmin=1e-2, vmax=1e2))
    cbar = plt.colorbar(strm, cax = cax, extend='both', orientation='horizontal')
    cbar.ax.yaxis.set_tick_params(pad=20)
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.tick_params(labelsize=global_fontsize)
    ax2.tick_params(axis='both', labelsize=global_fontsize)
    ax2.set_aspect(1)
    ax2.set_xlim(xmin, xmax)
    ax2.set_ylim(ymin, ymax)
    ax2.set_xlabel(r'$x$, [$c/\omega_{pl}$]', fontsize=global_fontsize)
#     ax2.set_yticks([])
    props = dict(boxstyle='square', facecolor='white', alpha=0.9, edgecolor='none')
    ax2.text(0.05, 0.97, r'$B/B_{\rm up}$', transform=ax2.transAxes, fontsize=global_fontsize, verticalalignment='top', bbox=props)
    ax2.set_xticks([-4000,0,4000])

    min_e = 0.02
    max_e = 10.
    min_n = 1e3
    max_n = 1e8
    cnts, bns = np.histogram(plasma.g / 10000., bins=np.logspace(np.log10(min_e), np.log10(max_e), 300))
    cnts = cnts * stride
    bns = average(bns)
    # cnts *= bns # for: gamma d f(gamma) / d gamma

    ax3.plot(bns, cnts, color = 'blue', label = 'plasma', linewidth = 0.8)
    ax3.set_xscale('log')
    ax3.set_yscale('log', nonposy='clip')
    ax3.yaxis.tick_left()
    ax3.yaxis.set_label_position("left")

    ax3.legend(loc='upper center', fontsize=global_fontsize)
    ax3.ticklabel_format(fontsize=global_fontsize)

    ax3.set_xlim(min_e, max_e)
    ax3.set_ylim(min_n, max_n)
    ax3.set_xlabel(r'$\gamma / \sigma$, $[m_e c^2]$', fontsize=global_fontsize)
    ax3.set_ylabel(r'$\gamma~\mathrm{d}f(\gamma)/\mathrm{d}\gamma$', fontsize=global_fontsize)
    ax3.tick_params(axis='both', labelsize=global_fontsize)

#     ax3 = plot_spectrum(ax3, plasma.g / 10000.,
#                         stride=stride,
#                         label='plasma', color = 'blue', min_e = 0.02, max_e = 10, min_n = 1e3, max_n = 1e8)
#     ax3.set_xlabel(r'$\varepsilon / \sigma$, $[m_ec^2]$')

    plt.tight_layout()
    plt.savefig(output_dir + "all_" + str(step).zfill(3) + ".png", dpi=150)
