# plot cool_IC

import numpy as np

from aux_11 import *
from parser import define_variables


root = '/u/hhakoby1/outputs/cool_ic/'
output_dir = root + 'pics/'

simulation_variables = define_variables(root)

import os
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

code_downsampling = simulation_variables['code_downsampling']
skin_depth = 0.25
ppc0 = simulation_variables['ppc0']
speed_of_light = simulation_variables['speed_of_light']
stride = simulation_variables['stride']

step


plasma = getPlasma(root, step)

yglob_mid = (plasma.y.max() + plasma.y.min()) * 0.5
xglob_mid = (plasma.x.max() + plasma.x.min()) * 0.5

min_e = 0.02
max_e = 10.
min_n = 1e2
max_n = 1e7
cnts, bns = np.histogram(plasma.g / 10000., bins=np.logspace(np.log10(min_e), np.log10(max_e), 300))
cnts = cnts * stride
bns = average(bns)

np.savetxt(root + 'spec_{}.out'.format(step), (bns, cnts))
