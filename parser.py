import re

def define_variables(root):
    filename = root + 'input'
    tau0 = None
    bw_dens_lim = None
    cool_dens_lim = None
    epsph_max = None
    epsph_min = None
    lorrad_lap = None
    gamma_rad = None
    gamma_c = None
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('c '):
                speed_of_light = float(re.findall('\d+\.\d*', line)[0])
            if line.startswith('interval '):
                output_period = float(re.findall('\d+', line)[0])
            if line.startswith('istep '):
                code_downsampling = float(re.findall('\d+', line)[0])
            if line.startswith('stride '):
                stride = float(re.findall('\d+', line)[0])
            if line.startswith('sigma '):
                sigma = float(re.findall('\d+\.\d*', line)[0])
            if line.startswith('ppc0 '):
                ppc0 = float(re.findall('\d+', line)[1])
            if line.startswith('c_omp '):
                skin_depth = float(re.findall('\d+', line)[0])
            if line.startswith('bc_transfer_lap'):
                bc_transfer_lap = float(re.findall('\d+', line)[0])
            if line.startswith('lorc '):
                gamma_c = float(re.findall('\d+\.\d*', line)[0])
            if line.startswith('lorrad '):
                gamma_rad = float(re.findall('\d+\.\d*', line)[0])
            if line.startswith('tau0 '):
                tau0 = float(re.findall('\d+\.\d*', line)[0])
            if re.findall('lorrad_lap\s', line):
                lorrad_lap = float(re.findall('\d+', line)[0])
            if re.findall('epsph_min\s', line):
                epsph_min = float(re.findall('\d+\.?\d*', line)[0])
            if re.findall('epsph_max\s', line):
                epsph_max = float(re.findall('\d+\.?\d*', line)[0])
            if re.findall('bw_dens_lim\s', line):
                bw_dens_lim = float(re.findall('\d+', line)[0])
            if re.findall('cool_dens_lim\s', line):
                cool_dens_lim = float(re.findall('\d+', line)[0])
    return {
        'speed_of_light': speed_of_light,
        'output_period': output_period,
        'code_downsampling': code_downsampling,
        'stride': stride,
        'sigma': sigma,
        'ppc0': ppc0,
        'skin_depth': skin_depth,
        'bc_transfer_lap': bc_transfer_lap,
        'gamma_c': gamma_c,
        'gamma_rad': gamma_rad,
        'tau0': tau0,
        'lorrad_lap': lorrad_lap,
        'epsph_max': epsph_max,
        'epsph_min': epsph_min,
        'bw_dens_lim': bw_dens_lim,
        'cool_dens_lim': cool_dens_lim
    }
