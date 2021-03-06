import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import copy
import matplotlib.ticker as ticker
import scipy as sp
import scipy.interpolate
from matplotlib.ticker import NullFormatter
import os

from aux_11 import *
from color_data import plasma_cmap
from color_data import viridis_cmap
from color_data import magma_cmap
from color_data import inferno_cmap
plt.register_cmap(name='plasma', cmap=plasma_cmap)
plt.register_cmap(name='viridis', cmap=viridis_cmap)
plt.register_cmap(name='magma', cmap=magma_cmap)
plt.register_cmap(name='inferno', cmap=inferno_cmap)

def plot_dens(ax, x, y,
                  dens, vmin, vmax,
                  label,
                  xmin, xmax, ymin, ymax,
                  cmap, scaling, setover = None, setunder = None, extend = 'neither', ret_cbar = False):
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
        strm = ax.imshow(dens, cmap=my_cmap,
                         norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax),
                         extent=[x.min(), x.max(), y.min(), y.max()])
    elif scaling == 'lin':
        strm = ax.imshow(dens, cmap=my_cmap,
                         vmin=vmin, vmax=vmax,
                         extent=[x.min(), x.max(), y.min(), y.max()])
    elif scaling == 'symlog':
        strm = ax.imshow(dens, cmap=my_cmap,
                         norm=mpl.colors.SymLogNorm(linthresh=vmax/10., vmin=vmin, vmax=vmax),
                         extent=[x.min(), x.max(), y.min(), y.max()])
    cbar = plt.colorbar(strm, cax = cax, extend=extend)
    cbar.ax.yaxis.set_tick_params(pad=10)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    # cbar.ax.tick_params(labelsize=fontsize)
    ax.tick_params(axis='both')
    ax.set_aspect(1)
    ax.set_ylabel(r'$x$, [$c/\omega_{pl}$]')
    props = dict(boxstyle='square', facecolor='white', alpha=0.9, edgecolor='none')
    ax.text(0.02, 0.95, label, transform=ax.transAxes, verticalalignment='top', bbox=props)
    if ret_cbar:
        return (cbar, ax)
    else:
        return ax

def rainbow_fill(ax,X,Y,cmap='jet', alpha=1.):
    def rect(x,y,w,h,c):
        c = list((np.array(c) * alpha + 1. * (1. - alpha))[:-1])
        polygon = plt.Rectangle((x,y),w,h,color=c)
        ax.add_patch(polygon)
    dxs = (X[1:]-X[:-1])
    N  = float(X.size)
    for n,(x,y,dx) in enumerate(zip(X,Y,dxs)):
        color = plt.get_cmap(cmap)(n/N)
        rect(x-dx,0,dx,y,color)
    return ax

def plot_spectrum(ax, bins, cnts, nprocs, bin_size = 151,
                      label = None, color = 'black', ls = '-', lw = 0.5, ncol = 3,
                      min_e = 1e-1, max_e = 1e3, min_n = 1e0, max_n = 1e10, normalize = False, return_data = False):

    bins = bins[:bin_size]
    bins = np.append([1e-1], bins)[:-1]
    # bins = average_array(bins)

    cnts = reduce_array(cnts, nprocs, bin_size)

    indices = (bins > min_e) & (bins < max_e)
    bins = bins[indices]
    cnts = cnts[indices]

    if max(cnts) < min_n:
        cnts += min_n / 10.
    if normalize:
        cnts /= max(cnts)
    ax.plot(bins, cnts, c=color, ls=ls, label=label, lw=lw)
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.yaxis.tick_left()
    ax.yaxis.set_label_position("left")

    ax.legend(loc='upper center', ncol=ncol)
    # ax.ticklabel_format(fontsize=fontsize)

    ax.set_xlim(min_e, max_e)
    ax.set_ylim(min_n, max_n)
    ax.set_xlabel(r'$\varepsilon$')
    ax.set_ylabel(r'$\varepsilon~\mathrm{d}f/\mathrm{d}\varepsilon$')
    ax.tick_params(axis='both')
    if return_data:
        return (ax, bins, cnts)
    else:
        return ax

def plot_temperature(ax, plasma,
                     xmin, xmax, ymin, ymax,
                     max_g, skin_depth = 10, dwn = 8):
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
              epsph_min, epsph_max):
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

def plot_photonB_vs_gamma(ax, root, step, sigma, gamma_c, emin=1e-2):
    fname = root + 'phst.tot.' + str(step).zfill(3)
    if not os.path.isfile(fname):
        data_gams = np.array([1e1])
        data_Bs = np.array([1e-3])
    else:
        data = h5py.File(fname, 'r')
        data_gams = data['gam'].value
        data_Bs = data['B'].value
    my_cmap = copy.copy(mpl.cm.get_cmap('jet'))
    my_cmap.set_bad(my_cmap(0))
    my_cmap.set_under(my_cmap(0))

    b_ax = np.linspace(-3,1,100)
    g_ax = np.linspace(0,4,100)
    if len(data_gams) > 10000:
        cnt = ax.hist2d(np.log10(data_gams), np.log10(data_Bs), bins=(g_ax, b_ax), norm=mpl.colors.LogNorm(), cmap=my_cmap);
    else:
        data_gams = np.array([1e1])
        data_Bs = np.array([1e-3])
        cnt = ax.hist2d(np.log10(data_gams), np.log10(data_Bs), bins=(g_ax, b_ax), norm=mpl.colors.LogNorm(vmin=1, vmax=1e5), cmap=my_cmap);
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    cbar = plt.colorbar(cnt[-1], cax = cax)
    cbar.set_label(r'\# of photons')

    g_ax, b_ax = np.meshgrid(g_ax, b_ax)
    epsph = (10**(g_ax) / gamma_c)**2 * 10**(b_ax) * np.sqrt(sigma / 1000.)
    levs = np.logspace(np.log10(emin), 3, (3-np.log10(emin)) + 1)
    clab = ax.contour(g_ax, b_ax, epsph, levels=levs, norm=mpl.colors.LogNorm(), colors='black');
    if len(levs > 5):
        clab = ax.contour(g_ax, b_ax, epsph, levels=levs[0::2], norm=mpl.colors.LogNorm(), colors='black');

    g_ticks = [0, 1, 2, 3, 4]
    ax.set_xticks(g_ticks)
    ax.set_xticklabels([r'$10^{{{0}}}$'.format(str(pow)) for pow in g_ticks])
    b_ticks = [-3, -2, -1, 0, 1]
    ax.set_yticks(b_ticks)
    ax.set_yticklabels([r'$10^{{{0}}}$'.format(str(pow)) for pow in b_ticks])

    fmt = ticker.LogFormatterMathtext()
    fmt.create_dummy_axis()

    clabs = ax.clabel(clab, fmt=fmt);

    [txt.set_bbox(dict(facecolor='white', edgecolor='black', pad=6)) for txt in clabs]

    ax.set_xlabel(r'particle $\gamma$')
    ax.set_ylabel(r'$B / B_{\rm up}$');
    ax.set_xlim(np.log10(5),np.log10(1e4))
    ax.set_ylim(-3,1)

    mpl.rcParams['hatch.color'] = (0,0,0,.2)
    xs = np.linspace(-1,4,5)
    ys = np.log10(emin * (1e3 / sigma)**0.5 * (gamma_c / 10**xs)**2)
    ax.fill_between(xs, -5, ys, hatch="//", linewidth=0.0, alpha=1.0, color='white')
    ax.fill_between(xs, -5, ys, hatch="//", linewidth=0.0, alpha=0.0)
    # txt = ax.text(1.3, -1.8, "not tracked",
    #                      color='black', horizontalalignment='center', verticalalignment='center', rotation=-45)
    # txt.set_bbox(dict(facecolor='white', alpha=1, edgecolor='none'));

def plot_e1_vs_e2(ax, root, step, emin = 0.01):
    nullfmt = NullFormatter()
    fname = root + 'prst.tot.' + str(step).zfill(3)
    if not os.path.isfile(fname):
        E1s = np.array([1e10])
        E2s = np.array([1e10])
        cosphis = np.array([1])
    else:
        data = h5py.File(fname, 'r')
        E1s = data['E1'].value
        E2s = data['E2'].value
        cosphis = data['cosph'].value
        E2s = E2s[(E1s != -2) & (E1s != 0)]
        cosphis = cosphis[(E1s != -2) & (E1s != 0)]
        E1s = E1s[(E1s != -2) & (E1s != 0)]
    if len(E1s) < 1e4:
        E1s = np.array([1e10])
        E2s = np.array([1e10])
        cosphis = np.array([1])
    indices = np.random.choice(np.arange(len(E1s)), 100000)
    E1s = E1s[indices]
    E2s = E2s[indices]
    cosphis = cosphis[indices]

    scat = ax.scatter(E1s, E2s, s=1, c=np.arccos(cosphis)*180/np.pi, edgecolor='none')
    ax.plot(np.logspace(-10,10,5), 1. / np.logspace(-10,10,5), c='red', ls='--')
    # ax.plot(np.logspace(-2,10,5), [1e-2] * 5, c='black', ls='--')
    # ax.plot([1e-2] * 5, np.logspace(-2,10,5), c='black', ls='--')
    cbaxes = inset_axes(ax, width="30%", height="3%", loc=1, borderpad=2)
    cbar = plt.colorbar(scat, cax=cbaxes, ticks=np.linspace(5, 178, 3), orientation='horizontal')
    cbar.ax.set_xticks([5, 91.5, 178])
    cbar.ax.set_xticklabels(['0', '90', '180'])
    cbaxes.xaxis.set_ticks_position('bottom')
    cbar.set_label(r'$\phi$')
    cbaxes.xaxis.set_label_position('top')
    ax.set_yscale('log', nonposy='clip')
    ax.set_xscale('log', nonposx='clip')
    ax.set_xlim((1e-3, 1e4))
    ax.set_ylim((1e-3, 1e4))
    ax.set_xlabel(r'$\epsilon_1/m_ec^2$')
    ax.set_ylabel(r'$\epsilon_2/m_ec^2$')
    if emin > 1e-3:
        ax.axvline(emin, c='black', ls='--')
        ax.axhline(emin, c='black', ls='--')
        mpl.rcParams['hatch.color'] = (0,0,0,.2)
        ax.fill_between([1e-5,emin,emin,1e5],[1e5,1e5,emin,emin], hatch="\\\\", linewidth=0.0, alpha=0.0)
        txt = ax.text(1, emin / 3.333, r'not tracked', color='black', horizontalalignment='center', verticalalignment='center')
        txt.set_bbox(dict(facecolor='white', alpha=1, edgecolor='none'))
    mpl.rcParams['hatch.color'] = (1,0,0,.2)
    ax.fill_between(np.logspace(-10,10,2), 1. / np.logspace(-10,10,2), hatch="//", linewidth=0.0, alpha=0.0)

def plot_energies(ax, stride, data, step):
    colors = ['#c1493c', '#426cce', '#51dbb1', '#dbd259']
    ax.plot([], [], color = colors[0], label = 'B-field')
    ax.plot([], [], color = colors[1], label = 'initial plasma')
    ax.plot([], [], color = colors[2], label = 'secondary pairs')
    ax.plot([], [], color = colors[3], label = 'photons')
    ax.stackplot(data[:,0]*500, data[:,1], data[:,2], data[:,3], data[:,4], colors=colors)
    leg = ax.legend(loc='lower left')
    leg.get_frame().set_facecolor('#FFFFFF')
    leg.get_frame().set_linewidth(0)
    leg.get_frame().set_alpha(1)
    ax.set_xlim(0, np.max(data[:,0]*500))
    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'energy budgets')
    ax.axvline(step*500, c='blue', ls='--')
