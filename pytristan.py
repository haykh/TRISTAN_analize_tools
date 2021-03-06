import numpy as np
import h5py
import os, os.path
import pandas as pd
from parser import define_variables
from tqdm import tqdm

def getPlasma(root, i):
    data2 = h5py.File(root + 'prtl.tot.' + str(i).zfill(3),'r')
    gammae = np.sqrt(1. + data2['ue'].value**2 + data2['ve'].value**2 + data2['we'].value**2)
    gammai = np.sqrt(1. + data2['ui'].value**2 + data2['vi'].value**2 + data2['wi'].value**2)
    inde = np.sign(data2['inde'].value) * (data2['inde'].value + data2['proce'].value * 1e9)
    indi = np.sign(data2['indi'].value) * (data2['indi'].value + data2['proci'].value * 1e9)
    plasma = pd.DataFrame({'x': np.concatenate((data2['xe'].value, data2['xi'].value)),
                           'y': np.concatenate((data2['ye'].value, data2['yi'].value)),
                           'u': np.concatenate((data2['ue'].value, data2['ui'].value)),
                           'v': np.concatenate((data2['ve'].value, data2['vi'].value)),
                           'w': np.concatenate((data2['we'].value, data2['wi'].value)),
                           'g': np.concatenate((gammae, gammai)),
                           'ind': np.concatenate((inde, indi)),
                           'ch': np.concatenate(([-1] * len(data2['xe'].value), [1] * len(data2['xi'].value)))
                          })
    return plasma

def getPhotons(root, i):
    if os.path.isfile(root + 'phot.tot.' + str(i).zfill(3)):
        data1 = h5py.File(root + 'phot.tot.' + str(i).zfill(3),'r')
        energy = np.sqrt(data1['up'].value**2 + data1['vp'].value**2 + data1['wp'].value**2)
        phots = pd.DataFrame({'ch': data1['chp'].value,
                              'x': data1['xp'].value,
                             'y': data1['yp'].value,
                             'u': divideArray(data1['up'].value, energy),
                             'v': divideArray(data1['vp'].value, energy),
                             'w': divideArray(data1['wp'].value, energy),
                             'e': energy,
                             'ind': data1['indp'].value})
        phots = phots[phots.e > 0]
    else:
        phots = None
    return phots

def getSizes(root, i):
    f1=h5py.File(root + 'flds.tot.' + str(i).zfill(3),'r')
    hh=f1['densi'].value
    hh1=hh[0,:,:]
    hh1=np.transpose(hh1)
    hh1=np.rot90(hh1)
    return {
        'sizex': len(hh1[0]),
        'sizey': len(hh1)
    }

def getField(root, i, param, sizes = None,
             xmin = 0, xmax = -1, ymin = 0, ymax = -1):
    if not sizes:
        sizes = getSizes(root, i)
    f1 = h5py.File(root + 'flds.tot.' + str(i).zfill(3),'r')
    if param in f1.keys():
        hh = f1[param].value
    else:
        hh = np.array(f1['dens']) * 0.0
    hh1=hh[0,:,:]
    hh1=np.transpose(hh1)
    hh1=np.rot90(hh1)
    hh1 = hh1[ymin:ymax,xmin:xmax]
    hh1 = np.rot90(hh1, 3)
    return hh1

def getNumberOfFiles(directory):
    x = 0
    for root, dirs, files in os.walk(directory):
        x += sum(['flds' in fn for fn in files])
    return x

def average_array(array):
	return (array[1:] + array[:-1]) * 0.5
def reduce_array(arr, nprocs, bin_size):
    return np.sum(np.reshape(arr, (nprocs, bin_size)), axis=0)
# def reshape_array(arr, bin_size):
#     return np.array([arr[i:i + bin_size] for i in range(0, len(arr), bin_size)][0])


def divideArray( a, b ):
    """ ignore / 0, div0( [-1, 0, 1], 0 ) -> [0, 0, 0] """
    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide( a, b )
        c[ ~ np.isfinite( c )] = 0  # -inf inf NaN
    return c

def trackEnergy(root, finstep):
    simulation_variables = define_variables(root)

    code_downsampling = simulation_variables['code_downsampling']
    skin_depth = simulation_variables['skin_depth']
    ppc0 = simulation_variables['ppc0']
    speed_of_light = simulation_variables['speed_of_light']
    output_period = simulation_variables['output_period']
    stride = simulation_variables['stride']
    root += 'output/'

    m_el = (speed_of_light / skin_depth)**2 * (1. / ppc0)
    magnetic_energy = []
    particle_energy = []
    pairs_energy = []
    photons_energy = []
    npart = []
    realnpart = []
    for step in range(0, finstep):
        print ((int)(100. * step / finstep), '%')
        # B-field
        data = h5py.File(root + 'flds.tot.{}'.format(str(step).zfill(3)), 'r')
        bx = data['bx'].value
        by = data['by'].value
        bz = data['bz'].value
        ex = data['ex'].value
        ey = data['ey'].value
        ez = data['ez'].value
        b_en = np.sum(bx**2 + by**2 + bz**2 + ex**2 + ey**3 + ez**2) * code_downsampling**2 / (2. * m_el * speed_of_light**2)

        # particles
        parts = getPlasma(root, step)
        data = h5py.File(root + 'prtl.tot.{}'.format(str(step).zfill(3)), 'r')
        gammae = (data['gammae'].value)[data['inde'].value > 0]
        gammae_prs = (data['gammae'].value)[data['inde'].value < 0]
        prtl_en = 2. * np.sum(gammae - 1.) * stride
        pair_en = 2. * np.sum(gammae_prs - 1.) * stride

        # photons
        data = h5py.File(root + 'phot.tot.{}'.format(str(step).zfill(3)), 'r')
        u = data['up'].value
        v = data['vp'].value
        w = data['wp'].value
        ch = data['chp'].value
        phot_en = np.sqrt(u**2 + v**2 + w**2) * ch
        phot_en = np.sum(phot_en) * stride

        npart.append(len(data['up'].value) * stride)
        realnpart.append(np.sum(data['chp'].value) * stride)

        magnetic_energy.append(b_en)
        particle_energy.append(prtl_en)
        pairs_energy.append(pair_en)
        photons_energy.append(phot_en)
    return (np.array(magnetic_energy),
            np.array(particle_energy),
            np.array(pairs_energy),
            np.array(photons_energy),
            np.array(npart),
            np.array(realnpart))

def trackTimestep(root):
    filename = root + 'report'
    total = []
    laps = []
    with open(filename, 'r') as f:
        for line in f:
            if 'lap ' in line:
                laps.append((float)(re.findall("\d+", line)[0]))
            if 'Total, sec' in line:
                temp = (float)(re.findall("\d\.\d\d..\d\d", line)[0])
                total.append(temp)
    return (np.array(laps), np.array(total))

def track_energy(root, step, x_lim=200):
    simulation_variables = define_variables(root)
    code_downsampling = simulation_variables['code_downsampling']
    skin_depth = simulation_variables['skin_depth']
    ppc0 = simulation_variables['ppc0']
    speed_of_light = simulation_variables['speed_of_light']
    output_period = simulation_variables['output_period']
    stride = simulation_variables['stride']

    root += 'output/'
    sim_vals = h5py.File(root + 'param.000','r')
    mx0 = sim_vals['mx0'].value[0]
    my0 = sim_vals['my0'].value[0]

    def x_to_sd(x):
        return (x - mx0 * 0.5) / skin_depth
    def y_to_sd(y):
        return (y - my0 * 0.5) / skin_depth

    m_el = (speed_of_light / skin_depth)**2 * (1. / ppc0)

    # B-field
    bx = getField(root, step, 'bx', getSizes(root, step), ymin = 0, ymax = -1)
    by = getField(root, step, 'by', getSizes(root, step), ymin = 0, ymax = -1)
    bz = getField(root, step, 'bz', getSizes(root, step), ymin = 0, ymax = -1)
    ex = getField(root, step, 'ex', getSizes(root, step), ymin = 0, ymax = -1)
    ey = getField(root, step, 'ey', getSizes(root, step), ymin = 0, ymax = -1)
    ez = getField(root, step, 'ez', getSizes(root, step), ymin = 0, ymax = -1)
    x = (np.arange(len(bx[0])) - max(np.arange(len(bx[0]))) * 0.5) * code_downsampling / skin_depth
    y = (np.arange(len(bx)) - max(np.arange(len(bx))) * 0.5) * code_downsampling / skin_depth
    x, y = np.meshgrid(x, y)
    bx = np.where(np.abs(y)<x_lim, bx, 0)
    by = np.where(np.abs(y)<x_lim, by, 0)
    bz = np.where(np.abs(y)<x_lim, bz, 0)
    ex = np.where(np.abs(y)<x_lim, ex, 0)
    ey = np.where(np.abs(y)<x_lim, ey, 0)
    ez = np.where(np.abs(y)<x_lim, ez, 0)
    b_en = np.sum(bx**2 + by**2 + bz**2 + ex**2 + ey**3 + ez**2) * code_downsampling**2 / (2. * m_el * speed_of_light**2)

    # particles
    # data = h5py.File(root + 'prtl.tot.{}'.format(str(step).zfill(3)), 'r')
    prtls = getPlasma(root, step)
    prtls = prtls[x_to_sd(prtls.x) < x_lim]
    prs = prtls[prtls.ind < 0]
    prtls = prtls[prtls.ind > 0]
    prtl_en = np.sum(prtls.g) * stride
    prs_en = np.sum(prs.g) * stride

    # photons
    import os.path
    if os.path.isfile(root + 'phot.tot.{}'.format(str(step).zfill(3))):
        data = h5py.File(root + 'phot.tot.{}'.format(str(step).zfill(3)), 'r')
        u = data['up'].value
        v = data['vp'].value
        w = data['wp'].value
        ch = data['chp'].value
        ee = np.sqrt(u**2 + v**2 + w**2)
        indexes = (ee > 1)
        ee = ee[indexes]
        ch = ch[indexes]
        phot_en = np.sum(ee * ch) * stride
    else:
        phot_en = 0

    return (b_en,
            prtl_en, prs_en,
            phot_en)

def get_energies_vs_time(root, x_lim, steps, norm=True):
    data = []
    for step in tqdm(steps):
        b_en, prtl_en, prs_en, phot_en = track_energy(root, step, x_lim)
        if norm:
            tot_en = b_en + prtl_en + prs_en + phot_en
        else:
            tot_en = 1
        data.append([step, b_en/tot_en, prtl_en/tot_en, prs_en/tot_en, phot_en/tot_en])
    return np.array(data)


def getSpectrum(root, step, bin_size = 151, get_new_parts = False):
    parts_ = np.zeros(bin_size)
    phots_ = np.zeros(bin_size)
    new_parts = np.zeros(bin_size)

    param = h5py.File(root + 'param.{}'.format(str(step).zfill(3)), 'r')
    nprocs = param['sizey'].value[0] * param['sizex'].value[0]
    data = h5py.File(root + 'spec.tot.{}'.format(str(step).zfill(3)), 'r')
    bins = data['bn'].value
    parts = data['npart'].value
    pairs = data['npair'].value
    if pairs.max() > 0:
        parts = pairs
    phots = data['ninst'].value
    bins = bins[:bin_size-1]
    bins = np.append(10**np.floor(np.log10(bins.min())), bins)
    parts = reduce_array(parts, nprocs, bin_size)
    phots = reduce_array(phots, nprocs, bin_size)

    parts_ = np.array(parts)
    phots_ = np.array(phots)

    param = h5py.File(root + 'param.{}'.format(str(step).zfill(3)), 'r')
    nprocs = param['sizey'].value[0] * param['sizex'].value[0]

    prstfname = root + 'prst.tot.{}'.format(str(step).zfill(3))
    if os.path.isfile(prstfname) and get_new_parts:
        data = h5py.File(prstfname, 'r')
        if len(data['x'].value) > nprocs:
            data = data['x'].value
            gams = np.sqrt(data[data > 1])
            new_parts, bns = np.histogram(gams, bins=bins)
            new_parts = np.append(new_parts, 0)
    return {'bins': bins,
            'parts': parts_,
            'phots': phots_,
            'new_parts': np.array(new_parts)}

def getSpec(root, step, keywords = ('bn', 'npart', 'nphot')):
    param = h5py.File(root + 'param.{}'.format(str(step).zfill(3)), 'r')
    nprocs = param['sizey'].value[0] * param['sizex'].value[0]
    data = h5py.File(root + 'spec.tot.{}'.format(str(step).zfill(3)), 'r')
    bin_size = int(len(data['bn'].value) / nprocs)
    datas = []
    for key in keywords:
        data_temp = data[key].value
        if key == keywords[0]:
            data_temp = data_temp[:bin_size-1]
            data_temp = np.append(10**np.floor(np.log10(data_temp.min())), data_temp)
        else:
            data_temp = reduce_array(data_temp, nprocs, bin_size)
        datas.append(data_temp)
    datas = np.array(datas)
    return {key: value for (key, value) in zip(keywords, datas)}

def scroll_images(root, field, steps, istep = 4,
                  cmap='plasma', vmin=0.1, vmax=200,
                  overplot = lambda ax: ax,
                  draw_procs = False
                  ):
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import numpy as np
    from copy import copy
    import h5py
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib as mpl
    cmap = copy(plt.get_cmap(cmap))
    cmap.set_bad(cmap(0))
    n = len(steps)
    flds = []
    for i in range(n):
        step = steps[i]
        data = h5py.File(root + 'flds.tot.{}'.format(str(step).zfill(3)), 'r')
        flds.append(data[field].value[0])

    def remove_keymap_conflicts(new_keys_set):
        for prop in plt.rcParams:
            if prop.startswith('keymap.'):
                keys = plt.rcParams[prop]
                remove_list = set(keys) & new_keys_set
                for key in remove_list:
                    keys.remove(key)

    def multi_slice_viewer(volume):
        remove_keymap_conflicts({'j', 'k'})
        fig, ax = plt.subplots()
        ax.volume = volume
        ax.index = 0
        xmin = 0
        xmax = len(volume[0][0]) * istep
        ymin = 0
        ymax = len(volume[0]) * istep
        im = ax.imshow(volume[0],
                       norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax),
                       cmap=cmap, extent=(xmin,xmax,ymin,ymax))
        ax.set_title(steps[0])

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax);

        if draw_procs:
            cpus = h5py.File(root + 'cpus.box.{}'.format(str(steps[ax.index]).zfill(3)), 'r')
            ranks = cpus['rank'].value
            xmins = cpus['xmin'].value
            sizexs = cpus['xmax'].value - xmins
            ymins = cpus['ymin'].value
            sizeys = cpus['ymax'].value - ymins
            for rnk, xm, sx, ym, sy in zip(ranks, xmins, sizexs, ymins, sizeys):
                p = patches.Rectangle((xm, ym), sx, sy, linewidth=0.2, edgecolor='w', facecolor='none')
                ax.add_artist(p)

        ax = overplot(ax)

        fig.canvas.mpl_connect('key_press_event', process_key)

    def process_key(event):
        fig = event.canvas.figure
        ax = fig.axes[0]
        if event.key == 'j':
            previous_slice(ax)
        elif event.key == 'k':
            next_slice(ax)
        fig.canvas.draw()

    def previous_slice(ax):
        volume = ax.volume
        ax.index = (ax.index - 1) % n  # wrap around using %
        xmin = 0
        xmax = len(volume[ax.index][0]) * istep
        ymin = 0
        ymax = len(volume[ax.index]) * istep
        ax.images[0].set_array(volume[ax.index])
        ax.images[0].set_extent((xmin,xmax,ymin,ymax))
        ax.set_title(steps[ax.index])
        if draw_procs:
            cpus = h5py.File(root + 'cpus.box.{}'.format(str(steps[ax.index]).zfill(3)), 'r')
            ranks = cpus['rank'].value
            xmins = cpus['xmin'].value
            sizexs = cpus['xmax'].value - xmins
            ymins = cpus['ymin'].value
            sizeys = cpus['ymax'].value - ymins
            for rnk, xm, sx, ym, sy, art in zip(ranks, xmins, sizexs, ymins, sizeys, ax.artists):
                art.set_xy((xm, ym))
                art.set_height(sy)
                art.set_width(sx)

    def next_slice(ax):
        volume = ax.volume
        ax.index = (ax.index + 1) % n
        xmin = 0
        xmax = len(volume[ax.index][0]) * istep
        ymin = 0
        ymax = len(volume[ax.index]) * istep
        ax.images[0].set_array(volume[ax.index])
        ax.images[0].set_extent((xmin,xmax,ymin,ymax))
        ax.set_title(steps[ax.index])
        if draw_procs:
            cpus = h5py.File(root + 'cpus.box.{}'.format(str(steps[ax.index]).zfill(3)), 'r')
            ranks = cpus['rank'].value
            xmins = cpus['xmin'].value
            sizexs = cpus['xmax'].value - xmins
            ymins = cpus['ymin'].value
            sizeys = cpus['ymax'].value - ymins
            for rnk, xm, sx, ym, sy, art in zip(ranks, xmins, sizexs, ymins, sizeys, ax.artists):
                art.set_xy((xm, ym))
                art.set_height(sy)
                art.set_width(sx)

    multi_slice_viewer(flds)


def drawPanels(root, panels, steps,
         istep = 4,
         mykeys = ['j', 'k'], figsize=(5, 3)):
    import h5py
    import numpy as np
    from copy import copy
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib as mpl
    from matplotlib import rc
    rc('font', **{'family': 'serif', 'serif': ['Helvetica'], 'size': 15})
    rc('text', usetex=True)
    plt.style.use('fast')

    ncols = 2
    nrows = int(np.ceil(len(panels) / ncols))
    figsize=(figsize[0] * ncols, figsize[1] * nrows)

    mykeys_dict = {mykeys[0], mykeys[1]}
    n = len(steps)
    fields = dict()
    cmaps = dict()
    normfuncs = dict()
    vmins = dict()
    vmaxs = dict()

    def read_data():
        # read data for all steps and save to dict
        for step in steps:
            data = h5py.File(root + 'flds.tot.{}'.format(str(step).zfill(3)), 'r')
            for panel in panels:
                fld_name = panel['name']
                value = panel['value']
                if type(value) == type('string'):
                    new_data = data[value].value[0]
                else:
                    new_data = value(data)
                try:
                    fields[fld_name].append(new_data)
                except:
                    fields[fld_name] = [new_data]
                if step == steps[0]: # do only once
                    # read colormaps
                    try:
                        cmap = panel['cmap']
                    except:
                        cmap = 'viridis'
                    cmap = copy(plt.get_cmap(cmap))
                    cmap.set_bad(cmap(0))
                    cmaps[fld_name] = cmap
                    # read boundaries
                    try:
                        vmin = panel['vmin']
                    except:
                        vmin = 1e-1
                    vmins[fld_name] = vmin
                    try:
                        vmax = panel['vmax']
                    except:
                        vmax = 200.
                    vmaxs[fld_name] = vmax
                    # read norms
                    try:
                        norm = panel['norm']
                    except:
                        norm = 'log'
                    if norm == 'log':
                        normfunc = mpl.colors.LogNorm
                    else:
                        normfunc = mpl.colors.Normalize
                    normfuncs[fld_name] = normfunc

    def remove_keymap_conflicts(new_keys_set):
        for prop in plt.rcParams:
            if prop.startswith('keymap.'):
                keys = plt.rcParams[prop]
                remove_list = set(keys) & new_keys_set
                for key in remove_list:
                    keys.remove(key)

    def redraw(ax):
        volume = ax.volume
        xmin = 0
        xmax = len(volume[ax.index][0]) * istep
        ymin = 0
        ymax = len(volume[ax.index]) * istep
        ax.images[0].set_array(volume[ax.index])
        ax.images[0].set_extent((xmin,xmax,ymin,ymax))

    def initialize():
        remove_keymap_conflicts(mykeys_dict)
        fig, axes = plt.subplots(nrows, 2, figsize=figsize)
        for ii in range(ncols * nrows):
            if nrows == 1:
                ax = axes[ii]
            else:
                ax = axes[int(np.floor(ii / ncols))][ii % ncols]
            if ii >= len(panels):
                panel = panels[-1]
                fld_name = ''
                field = fields[panel['name']]
                im = ax.imshow(np.array(field[0]) * 0.)
                ax.update = False
            else:
                panel = panels[ii]
                fld_name = panel['name']
                panel = panels[ii]
                vmin = vmins[fld_name]
                vmax = vmaxs[fld_name]
                cmap = cmaps[fld_name]
                field = fields[fld_name]
                normfunc = normfuncs[fld_name]
                im = ax.imshow(field[0], norm=normfunc(vmin=vmin, vmax=vmax), cmap=cmap)
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                cbar = plt.colorbar(im, cax=cax, label=fld_name)
                ax.update = True
            ax.set_aspect(1)
            ax.volume = field
            ax.index = 0
            redraw(ax)
        fig.index = 0
        fig.suptitle(steps[fig.index])
        fig.canvas.mpl_connect('key_press_event', process_key)

    def process_key(event):
        fig = event.canvas.figure
        if event.key == mykeys[0]:
            fig.index = (fig.index - 1) % n
            fig.suptitle(steps[fig.index])
            for ax in fig.axes:
                if ax.update:
                    prev_step(ax)
        elif event.key == mykeys[1]:
            fig.index = (fig.index + 1) % n
            fig.suptitle(steps[fig.index])
            for ax in fig.axes:
                if ax.update:
                    next_step(ax)
        fig.canvas.draw()

    def prev_step(ax):
        ax.index = (ax.index - 1) % n
        redraw(ax)

    def next_step(ax):
        ax.index = (ax.index + 1) % n
        redraw(ax)

    read_data()
    initialize()


def reconnectionDiagnostics(root, steps, parameters=None, unit_convert=True, plt_size=(15,8)):
    import numpy as np
    from tqdm import tqdm_notebook
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    from scipy.optimize import curve_fit
    from matplotlib import rc

    try:
        L = parameters['L']
    except:
        L = 14000.
    try:
        sigma = parameters['sigma']
    except:
        sigma = 10.
    try:
        interval = parameters['interval']
    except:
        interval = 500
    try:
        sd = parameters['skin_depth']
    except:
        sd = 5.
    rL = np.sqrt(sigma) * sd
    istep = 4
    c = 0.45
    g_star = 20
    if unit_convert:
        step_to_time = interval / (rL / c)
        cell_to_rl = istep / rL
    else:
        step_to_time = 1
        cell_to_rl = istep

    rc('font', **{'family': 'serif', 'serif': ['Helvetica'], 'size': 10})
    rc('text', usetex=True)
    plt.style.use('fast')

    # density
    def get_profile(root, steps):
        flds = []
        for step in tqdm_notebook(steps, leave=False):
            data = h5py.File(root + 'flds.tot.{}'.format(str(step).zfill(3)))
            flds.append(np.mean(data['dens'].value[0], axis=1))
        return np.array(flds)

    data1 = get_profile(root, steps)

    Y = steps * step_to_time
    X = np.arange(len(data1[0])) * cell_to_rl

    # inflow
    def get_inflows(root, steps):
        infs = []
        for step in tqdm_notebook(steps, leave=False):
            data = h5py.File(root + 'flds.tot.{}'.format(str(step).zfill(3)))
            beta_x = data['v3x'].value[0]
            mid_ = int(len(beta_x[0]) / 2)
            beta_x = np.abs(beta_x[:, mid_ - int(mid_ * 0.2) : mid_ + int(mid_ * 0.2)])
            infs.append(np.mean(beta_x))
        return np.array(infs)
    data2 = get_inflows(root, steps)

    plt.figure(figsize=plt_size)
    ax_dens = plt.subplot2grid((3, 2), (0, 0), rowspan=2)
    ax_spec = plt.subplot(322)
    ax_gam = plt.subplot(324)
    ax_pow = plt.subplot(326)
    ax_inf = plt.subplot(325)

    ax_dens.pcolorfast(X, Y, data1, cmap='plasma');
    if unit_convert:
        ax_dens.set_xlabel(r'$x$ [$r_L$]')
        ax_dens.set_ylabel(r'$t$ [$r_{L}/c$]')
    else:
        ax_dens.set_xlabel(r'$x$ [cells]')
        ax_dens.set_ylabel(r'$t$ [steps]')
    ax_dens.set_xlim(X.min(), X.max())
    ax_dens.set_ylim(Y.min(), Y.max())
    ax2_dens = ax_dens.twinx()
    ax2_dens.set_ylabel(r"steps")
    ax2_dens.set_ylim(steps.min(), steps.max())

    ax_inf.plot(Y, data2, c='red')
    if unit_convert:
        ax_inf.set_xlabel(r'$ct/r_{\rm L}$')
    else:
        ax_inf.set_xlabel(r'$t$ [steps]')
    ax_inf.set_ylabel(r'$\beta_{\rm in}$')
    ax_inf.set_xscale('log')
    ax_inf.set_xlim(1e2, max(Y)*1.2)

    # spectrum
    if len(steps) > 50:
        inter = int(len(steps) / 50.)
        steps_spec = steps[0:-1:inter]
    else:
        steps_spec = steps
    times = steps_spec * interval / (rL / c)
    cols = plt.cm.viridis(np.linspace(0, 1, len(steps_spec)))
    for ii, step in zip(range(len(steps_spec)), steps_spec):
        spec = getSpec(root, step)
        im = ax_spec.plot(spec['bn'], spec['npart'], c=cols[ii])

    sm = plt.cm.ScalarMappable(cmap='viridis')
    sm.set_array(np.linspace(times.min(), times.max(), 10))
    cbaxes = inset_axes(ax_spec, width="40%", height="5%", loc='upper right',
                        bbox_to_anchor=(-0.01, 0.46, 1, 0.5), bbox_transform=ax_spec.transAxes)
    cb = plt.colorbar(sm, cax=cbaxes, orientation='horizontal')
    cb.ax.xaxis.set_ticks_position('bottom')
    cb.ax.xaxis.set_label_position('top')
    if unit_convert:
        cb.set_label(r'$ct/r_L$')
    else:
        cb.set_label(r'$t$ [steps]')


    ax_spec.set_ylim(1e2, 2e7)
    ax_spec.set_xlim(1e-1, 1e4)
    ax_spec.set_xscale('log')
    ax_spec.set_yscale('log');
    ax_spec.set_xlabel(r'$\gamma-1$')
    ax_spec.set_ylabel(r'$(\gamma-1)f(\gamma)$');

    # fits
    def func(x, n0, p, g_cut):
        return n0 * (x - 1) * (x / g_star)**(-p) * np.exp(-x / g_cut)
    tts = []
    gcuts = []
    gmaxs = []
    pows = []
    times = steps_spec * interval / (rL / c)
    for ii, step in zip(range(len(steps_spec)), steps_spec):
        spec = getSpec(root, step)
        xs = spec['bn']
        ys = spec['npart']
        indices = (xs + 1 > g_star)
        ys1 = ys[indices]
        xs1 = xs[indices]
        gmax = np.trapz((xs1 + 1)**5 * ys1 / xs1, xs1+1) / np.trapz((xs1 + 1)**4 * ys1 / xs1, xs1 + 1)
        gmaxs.append(gmax)
        popt, pcov = curve_fit(func, xs1, ys1, bounds=([0, 1, 2], [np.inf, 5, 1e3]))
        pows.append(popt[1])
        gcuts.append(popt[2])
        tts.append(times[ii])

    xs = tts[1:]
    ys = np.sqrt(tts[1:] / tts[3]) * (gcuts[3] + gmaxs[3])*0.5
    ax_gam.plot(xs, ys, c=[.7]*3, ls='--', zorder=-1)
    ys = (tts[1:] / tts[3]) * (gcuts[3] + gmaxs[3])*0.5
    ax_gam.plot(xs, ys, c=[.7]*3, ls=':', zorder=-1)

    ax_gam.scatter(tts, gcuts, label=r'$\gamma_{\rm cut}$', c='blue')
    ax_gam.scatter(tts, gmaxs, label=r'$\gamma_{\rm max}$', c='red')
    ax_gam.set_xscale('log')
    ax_gam.set_yscale('log')
    ax_gam.legend()
    if unit_convert:
        ax_gam.set_xlabel(r'$tc/r_L$')
    else:
        ax_gam.set_xlabel(r'$t$ [steps]')
    ax_gam.set_ylabel(r'$\gamma$');
    ax_gam.set_xlim(1e2, max(tts)*1.2)

    ax_pow.scatter(tts, pows, label=r'$p$', c='brown')
    xs = tts
    ys = np.array(tts)*0 + 2
    ax_pow.plot(xs, ys, c=[.7]*3, ls='--', zorder=-1)
    ax_pow.set_xlim(1e2, max(tts)*1.2)
    ax_pow.set_ylim(1.2, 2.5)
    if unit_convert:
        ax_pow.set_xlabel(r'$tc/r_L$')
    else:
        ax_pow.set_xlabel(r'$t$ [steps]')
    ax_pow.set_ylabel(r'$p$');
    ax_pow.set_xscale('log')

    plt.tight_layout();


# def determineMaxDensity(root, start, end, fld):
#     maximum = 0
#     sizes = getSizes(root, start)
#     for step in range(start, end):
#         dens = getField(root, step, fld, sizes)
#         dens_max = dens.max()
#         if dens_max > maximum:
#             maximum = dens_max
#     return maximum
#
# def rebin(a, shape):
#     sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
#     return a.reshape(sh).mean(-1).mean(1)
#
# def determineMaxMultiplicity(root, start, end):
#     maximum = 0
#     sizes = getSizes(root, start)
#     for step in range(start, end):
#         dens = getField(root, step, 'dens', sizes)
#         denbw = getField(root, step, 'denbw', sizes)
#         dens = dens - denbw
#         multiplicity = divideArray(denbw, (dens - denbw))
#         h, w = multiplicity.shape
#         dh = h - 10*(h/10)
#         dw = w - 10*(w/10)
#         multiplicity = multiplicity[dh:,dw:]
#         h, w = multiplicity.shape
#         multiplicity = rebin(multiplicity, (h/10, w/10))
#         mul_max = multiplicity.max()
#         if mul_max > maximum:
#             maximum = mul_max
#     return maximum
#
# def determineMaxBsquared(root, start, end):
#     maximum = 0
#     sizes = getSizes(root, start)
#     for step in range(start, end):
#         bx = getField(root, step, 'bx', sizes)
#         by = getField(root, step, 'by', sizes)
#         bz = getField(root, step, 'bz', sizes)
#         bsqr = bx**2 + by**2 + bz**2
#         bsqr_max = bsqr.max()
#         if bsqr_max > maximum:
#             maximum = bsqr_max
#     return maximum
#
# def determineMaxVector(root, start, end, vector):
#     maximum = 0
#     sizes = getSizes(root, start)
#     for step in range(start, end):
#         vec = getField(root, step, vector, sizes)
#         vec_max = np.abs(vec).max()
#         if vec_max > maximum:
#             maximum = vec_max
#     return maximum
#
# def fmt(x):
#     a, b = '{:.0e}'.format(x).split('e')
#     b = int(b)
#     return r'${} \times 10^{{{}}}$'.format(a, b)
