import h5py
import pandas as pd
import numpy as np
import re

from aux_11 import *
from parser import define_variables

def trackEnergy(root, finstep, photQ, pairQ):
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
    pairs_energy = []
    npart = []
    realnpart = []
    for step in range(-1, finstep):
        print 100. * step / finstep
        bx = getField(root, step, 'bx', getSizes(root, step))
        by = getField(root, step, 'by', getSizes(root, step))
        bz = getField(root, step, 'bz', getSizes(root, step))
        b_en = (bx**2 + by**2 + bz**2) / (2. * m_el * speed_of_light**2)
        plasma = getPlasma(root, step)
        if photQ:
            photons = getPhotons(root, step)
        if pairQ:
            pairs = plasma[plasma.ind < 0]
            plasma = plasma[plasma.ind > 0]
        b_tot = np.sum(b_en) * code_downsampling**2
        prtl_tot = np.sum(plasma.g) * stride
        if pairQ:
            pair_tot = np.sum(pairs.g) * stride
        if photQ:
            phot_tot = np.sum(photons.e * photons.ch) * stride
            npart.append(len(photons) * stride)
            realnpart.append(np.sum(photons.ch) * stride)
        magnetic_energy.append(b_tot)
        particle_energy.append(prtl_tot)
        if pairQ:
            pairs_energy.append(pair_tot)
        if photQ:
            photons_energy.append(phot_tot)
    if not pairQ:
        pairs_energy = np.array(magnetic_energy) * 0.
    if not photQ:
        phot_energy = np.array(magnetic_energy) * 0.
        npart = phot_energy
        realnpart = npart
    return (np.array(magnetic_energy),
            np.array(particle_energy),
            np.array(photons_energy),
            np.array(pairs_energy),
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
                temp_pcent = float(re.findall("(\d+\.\d*)%", line)[0])
                temp += temp * temp_pcent / 100.
                total.append(temp)
    return (np.array(laps), np.array(total))

directory = '/u/hhakoby1/outputs/new_merging/'

mag_e, part_e, phot_e, prs_e, nph, rl_nph = trackEnergy(directory, 88, True, True)
laps, total = trackTimestep(directory)

np.savetxt('/u/hhakoby1/vis/new_energies.out', (mag_e, part_e, phot_e, prs_e, nph, rl_nph))
np.savetxt('/u/hhakoby1/vis/new_timesteps.out', (laps, total))
