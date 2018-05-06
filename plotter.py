import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import copy
import scipy as sp
import scipy.interpolate

from aux_11 import *
from color_data import plasma_cmap
from color_data import viridis_cmap
from color_data import magma_cmap
from color_data import inferno_cmap
plt.register_cmap(name='plasma', cmap=plasma_cmap)
plt.register_cmap(name='viridis', cmap=viridis_cmap)
plt.register_cmap(name='magma', cmap=magma_cmap)
plt.register_cmap(name='inferno', cmap=inferno_cmap)

global_fontsize = 15

def plot_dens(ax, x, y,
              dens, vmin, vmax,
              label,
              xmin, xmax, ymin, ymax,
              cmap, scaling, setover = None, setunder = None, extend = 'neither', ret_cbar = False,
              fontsize=global_fontsize):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    my_cmap = copy.copy(mpl.cm.get_cmap(cmap))
    if setunder is None:
        my_cmap.set_bad(my_cmap(0))
        my_cmap.set_under(my_cmap(0))
    else:
        my_cmap.set_bad(setunder)
        my_cmap.set_under(setunder)
    if setover is None:
        my_cmap.set_over(my_cmap(255))
    else:
        my_cmap.set_over(setover)
    if scaling == 'log':
        strm = ax.pcolorfast(x, y, dens, cmap=my_cmap, norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax))
    elif scaling == 'lin':
        strm = ax.pcolorfast(x, y, dens, cmap=cmap, vmin=vmin, vmax=vmax)
    elif scaling == 'symlog':
        strm = ax.pcolorfast(x, y, dens, cmap=cmap, norm=mpl.colors.SymLogNorm(linthresh=vmax/10., vmin=vmin, vmax=vmax))
    cbar = plt.colorbar(strm, cax = cax, extend=extend)
    cbar.ax.yaxis.set_tick_params(pad=10)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    cbar.ax.tick_params(labelsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize)
    ax.set_aspect(1)
    ax.set_ylabel(r'$x$, [$c/\omega_{pl}$]', fontsize=fontsize)
    props = dict(boxstyle='square', facecolor='white', alpha=0.9, edgecolor='none')
    ax.text(0.02, 0.95, label, transform=ax.transAxes, fontsize=fontsize, verticalalignment='top', bbox=props)
    if ret_cbar:
        return (cbar, ax)
    else:
        return ax

def plot_spectrum(ax, prtls, stride = 1,
                  label = None, color = 'black', ls = '-',
                  weights = None, min_e = 1e-1, max_e = 1e3, min_n = 1e0, max_n = 1e10, interp = False,
                  fontsize=global_fontsize, normalize = False):

    cnts, bns = np.histogram(prtls, bins=np.logspace(np.log10(min_e), np.log10(max_e), 150), weights = weights)
    cnts = cnts * stride
    bns = average(bns)
    # g df(g) / dg

    if normalize:
        cnts = np.array(cnts)
        cnts /= cnts.max()

    if len(bns[cnts != 0]) < 5:
        return ax

    if interp:
        # interpolation
        cnts += 1
        bns_new = np.logspace(np.log10(bns[0]), np.log10(bns[-1]), 1000)
        spl = sp.interpolate.splrep(np.log10(bns), np.log10(cnts), s=np.log10(cnts).max() / 10.)
        cnts_new = 10.0**(sp.interpolate.splev(np.log10(bns_new), spl))
        ax.plot(bns_new, cnts_new, color = color, ls = ls, label = label, linewidth = 0.8)
    else:
        ax.step(bns, cnts, color = color, ls = ls, label = label, linewidth = 0.8)


    ax.set_xscale('log')
    ax.set_yscale('log', nonposy='clip')
    ax.yaxis.tick_left()
    ax.yaxis.set_label_position("left")

    ax.legend(loc='upper center', ncol=5, fontsize=fontsize)
    ax.ticklabel_format(fontsize=fontsize)

    ax.set_xlim(min_e, max_e)
    ax.set_ylim(min_n, max_n)
    ax.set_xlabel(r'$\varepsilon$, $[m_e c^2]$', fontsize=fontsize)
    ax.set_ylabel(r'$\varepsilon~\mathrm{d}f(\varepsilon)/\mathrm{d}\varepsilon$', fontsize=fontsize)

        # ax.plot([1e2,1e4], [1e9, 1e7], color='purple', ls='--')
        # ax.text(2e3, 5e8, r'$\propto\gamma^{-1}$', fontsize=1.2*fontsize)
    ax.tick_params(axis='both', labelsize=fontsize)
    return ax

def plot_spectrum_new(ax, bins, cnts, nprocs, bin_size = 151,
                      label = None, color = 'black', ls = '-', lw = 0.5, ncol = 3,
                      min_e = 1e-1, max_e = 1e3, min_n = 1e0, max_n = 1e10,
                      fontsize=global_fontsize, normalize = False):
    def reduce_array(arr):
        return np.sum(np.reshape(arr, (nprocs, bin_size)), axis=0)
    def reshape_arr(arr):
        return np.array([arr[i:i + bin_size] for i in xrange(0, len(arr), bin_size)][0])

    bins = reshape_arr(bins)
    bins = np.append([1e-1], bins)
    bins = average(bins)

    cnts = reduce_array(cnts)

    indices = (bins > min_e) & (bins < max_e)
    bins = bins[indices]
    cnts = cnts[indices]

    if max(cnts) < min_n:
        cnts += min_n / 10.
    if normalize:
        cnts /= max(cnts)
    ax.plot(bins, cnts, c=color, ls=ls, label=label, lw=lw, drawstyle='steps')
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.yaxis.tick_left()
    ax.yaxis.set_label_position("left")

    ax.legend(loc='upper center', ncol=ncol, fontsize=fontsize)
    ax.ticklabel_format(fontsize=fontsize)

    ax.set_xlim(min_e, max_e)
    ax.set_ylim(min_n, max_n)
    ax.set_xlabel(r'$\varepsilon$', fontsize=fontsize)
    ax.set_ylabel(r'$\varepsilon~\mathrm{d}f/\mathrm{d}\varepsilon$', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize)
    return ax

def plot_temperature(ax, plasma,
                     xmin, xmax, ymin, ymax,
                     max_g, skin_depth = 10, dwn = 8,
                     fontsize=global_fontsize):
    dx = plasma.x.max() - plasma.x.min()
    dy = plasma.y.max() - plasma.y.min()
    cnts1, xed, yed = np.histogram2d(plasma.y, plasma.x, bins=(int(dy / dwn),int(dx / dwn)))
    cnts2, xed, yed = np.histogram2d(plasma.y, plasma.x, bins=(int(dy / dwn),int(dx / dwn)), weights=plasma.g)
    cnts = divideArray(cnts2, cnts1)

    cnts = np.transpose(cnts)

    x = np.arange(len(cnts[0])) * dwn / skin_depth
    x -= (x.max() + x.min())*0.5
    y = np.arange(len(cnts)) * dwn / skin_depth
    y -= (y.max() + y.min())*0.5

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    pcol = ax.pcolorfast(x, y, cnts, vmin = 0, vmax = max_g, cmap = 'inferno')
    cbar = plt.colorbar(pcol, cax = cax, extend='max')
    cbar.ax.tick_params(labelsize=fontsize)
    ax.set_aspect(1)
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(xmin, xmax)
    ax.tick_params(axis='both', labelsize=fontsize)
    ax.set_ylabel(r'$x$, [$c/\omega_{pl}$]', fontsize=fontsize)
    props = dict(boxstyle='square', facecolor='white', alpha=0.9, edgecolor='none')
    ax.text(0.02, 0.95, r'Average $\gamma$', transform=ax.transAxes, fontsize=fontsize, verticalalignment='top', bbox=props)
    return ax

def plot_stat(ax, root, step,
              epsph_min, epsph_max,
              fontsize=global_fontsize):
    if not os.path.isfile(root + 'stat.tot.' + str(step+1).zfill(3)):
        return ax
    data = h5py.File(root + 'stat.tot.' + str(step+1).zfill(3),'r')
    E1s = data['E1'].value
    E2s = data['E2'].value
    cosphis = data['cosph'].value
    E2s = E2s[(E1s != -2) & (E1s != 0)]
    cosphis = cosphis[(E1s != -2) & (E1s != 0)]
    E1s = E1s[(E1s != -2) & (E1s != 0)]

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    p = ax.scatter(E1s, E2s, s=5, c=np.arccos(cosphis)*180/np.pi)
    ax.plot(np.logspace(-4,4,100),1/np.logspace(-4,4,100), c='black', ls='--')
    cbar = plt.colorbar(p, cax = cax)
    cbar.set_label('relative angle', rotation=90, fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize)
    ax.set_xlabel(r'$\varepsilon_1$, [$m_ec^2$]', fontsize=fontsize)
    ax.set_ylabel(r'$\varepsilon_2$, [$m_ec^2$]', fontsize=fontsize)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(epsph_min, epsph_max)
    ax.set_ylim(epsph_min, epsph_max)
    ax.set_aspect(1)
    return ax
