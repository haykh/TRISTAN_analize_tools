import numpy as np
import h5py
import os, os.path
import pandas as pd
from aux_11 import *
from parser import define_variables


def getPlasma(root, i):
    data2 = h5py.File(root + 'prtl.tot.' + str(i).zfill(3),'r')
    gammae = np.sqrt(1. + data2['ue'].value**2 + data2['ve'].value**2 + data2['we'].value**2)
    gammai = np.sqrt(1. + data2['ui'].value**2 + data2['vi'].value**2 + data2['wi'].value**2)
    plasma = pd.DataFrame({'x': np.concatenate((data2['xe'].value, data2['xi'].value)),
                           'y': np.concatenate((data2['ye'].value, data2['yi'].value)),
                           'u': np.concatenate((data2['ue'].value, data2['ui'].value)),
                           'v': np.concatenate((data2['ve'].value, data2['vi'].value)),
                           'w': np.concatenate((data2['we'].value, data2['wi'].value)),
                           'g': np.concatenate((gammae, gammai)),
                           'ind': np.concatenate((data2['inde'].value, data2['indi'].value)),
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

def getField(root, i, param, sizes,
             xmin = 0, xmax = -1, ymin = 0, ymax = -1):
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


def average(array):
	newarr = array[:-1]
	for i in range(len(array) - 1):
		newarr[i] = 0.5 * (array[i] + array[i + 1])
	return newarr

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
        print int(100. * step / finstep), '%'
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
                laps.append(float(re.findall("\d+", line)[0]))
            if 'Total, sec' in line:
                temp = float(re.findall("\d\.\d\d..\d\d", line)[0])
                total.append(temp)
    return (np.array(laps), np.array(total))

def track_energy(root, step):
    simulation_variables = define_variables(root)
    code_downsampling = simulation_variables['code_downsampling']
    skin_depth = simulation_variables['skin_depth']
    ppc0 = simulation_variables['ppc0']
    speed_of_light = simulation_variables['speed_of_light']
    output_period = simulation_variables['output_period']
    stride = simulation_variables['stride']

    m_el = (speed_of_light / skin_depth)**2 * (1. / ppc0)
    root += 'output/'

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
    data = h5py.File(root + 'prtl.tot.{}'.format(str(step).zfill(3)), 'r')
    prtls = getPlasma(root, step)
    prs = prtls[prtls.ind < 0]
    prtls = prtls[prtls.ind > 0]
    prtl_en = np.sum(prtls.g - 1.) * stride
    prs_en = np.sum(prtls.g - 1.) * stride

    # photons
    data = h5py.File(root + 'phot.tot.{}'.format(str(step).zfill(3)), 'r')
    u = data['up'].value
    v = data['vp'].value
    w = data['wp'].value
    ch = data['chp'].value
    phot_en = np.sum(np.sqrt(u**2 + v**2 + w**2) * ch) * stride

    return (b_en,
            prtl_en, prs_en,
            phot_en)

def get_energies_vs_time(root, step_min=0, step_max=1):
    data = []
    for step in range(step_min, step_max + 1):
        b_en, prtl_en, prs_en, phot_en = track_energy(root, step)
        data.append([[step, b_en, prtl_en, prs_en, phot_en]])
    return data

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
