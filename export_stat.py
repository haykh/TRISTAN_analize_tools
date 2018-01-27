import sys
sys.path.append('/u/hhakoby1/vis/tqdm/')
from tqdm import tqdm

import h5py
import pandas as pd
import numpy as np
import re

from aux_11 import *
from parser import define_variables

def trackEnergy(root, finstep):
    simulation_variables = define_variables(root)

    code_downsampling = simulation_variables['code_downsampling']
    skin_depth = simulation_variables['skin_depth']
    ppc0 = simulation_variables['ppc0']
    speed_of_light = simulation_variables['speed_of_light']
    output_period = simulation_variables['output_period']
    stride = simulation_variables['stride']

    m_el = (speed_of_light / skin_depth)**2 * (1. / ppc0)
    magnetic_energy = []
    particle_energy = []
    photons_energy = []
    npart = []
    realnpart = []
    for step in tqdm(range(0, finstep)):
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
        data = h5py.File(root + 'prtl.tot.{}'.format(str(step).zfill(3)), 'r')
        prtl_en = 2. * np.sum(data['gammae'].value - 1.) * stride

        # photons
        data = h5py.File(root + 'phot.tot.{}'.format(str(step).zfill(3)), 'r')
        u = data['up'].value
        v = data['vp'].value
        w = data['wp'].value
        ch = data['chp'].value
        phot_en = np.sqrt(u**2 + v**2 + w**2) * ch
        # phot_en = phot_en[~ np.isnan(phot_en)]
        phot_en = np.sum(phot_en) * stride

        npart.append(len(data['up'].value) * stride)
        realnpart.append(np.sum(data['chp'].value) * stride)

        magnetic_energy.append(b_en)
        particle_energy.append(prtl_en)
        photons_energy.append(phot_en)
    return (np.array(magnetic_energy),
            np.array(particle_energy),
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

directory = '/u/hhakoby1/outputs/no_merging/'

mag_e, part_e, phot_e, nph, rl_nph = trackEnergy(directory, 45)
laps, total = trackTimestep(directory)

np.savetxt('/u/hhakoby1/vis/no_energies.out', (mag_e, part_e, phot_e, nph, rl_nph))
np.savetxt('/u/hhakoby1/vis/no_timesteps.out', (laps, total))
