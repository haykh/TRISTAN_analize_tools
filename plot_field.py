import h5py
from matplotlib.mlab import bivariate_normal
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors

from aux_11 import *
from color_data import plasma_cmap
plt.register_cmap(name='plasma', cmap=plasma_cmap)
from parser import define_variables
import helper as hlp

from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

name = 'test_02'
root = '/tigress/PERSEUS/hakobyan/TRISTAN/oct_2017/' + name + '/'
simulation_variables = define_variables(root)
root += 'output/'

output_dir = 'pics/october/' + name + '/'

max_number = int(getNumberOfFiles(root))
print "There are overall approx " + str(max_number) + " files."

code_downsampling = simulation_variables['code_downsampling']
skin_depth = simulation_variables['skin_depth']
ppc0 = simulation_variables['ppc0']
speed_of_light = simulation_variables['speed_of_light']
output_period = simulation_variables['output_period']
deposit_downsampling = 1.0

# maximum_dens = determineMaxDensity(root, 0, max_number - 1, 'dens')
# maximum_dens /= ppc0
#
# maximum_denbw = determineMaxDensity(root, 0, max_number - 1, 'denbw')
# maximum_denbw /= ppc0
#
# maximum_dens -= maximum_denbw

maximum_densph = determineMaxDensity(root, 0, max_number - 1, 'densph')
maximum_dens = maximum_densph

# start = int(raw_input("Start from: "))
# end = int(raw_input("End at: "))

start = 0
end = max_number

hlp.printProgress(0, 1, prefix = 'Progress:', barLength = 50)
for step in range(start, min(max_number, end)):
	# dens = getField(root, step, 'dens', getSizes(root, step))
	densph = getField(root, step, 'densph', getSizes(root, step))
	dens = densph
	# denbw = getField(root, step, 'denbw', getSizes(root, step))
	# dens -= denbw

	data_energies = getEnergies(root, step)
	photons = data_energies['phot']
	photons_ch = data_energies['ph_ch']
	# pairs_e = data_energies['lecs_p']
	# pairs_i = data_energies['ions_p']
	# pairs = np.concatenate((pairs_e, pairs_i), axis=0)
	# particles_e = data_energies['lecs']
	# particles_i = data_energies['ions']
	# particles = np.concatenate((particles_e, particles_i), axis=0)

	x = (np.arange(len(dens[0])) - max(np.arange(len(dens[0]))) * 0.5) * code_downsampling / skin_depth
	y = (np.arange(len(dens)) - max(np.arange(len(dens))) * 0.5) * code_downsampling / skin_depth
	x,y = np.meshgrid(x,y)

	fig = plt.figure(figsize=(15, 18))

	ax0 = plt.subplot2grid((2,1),(0,0))
	# ax1 = plt.subplot2grid((4,3),(1,0),colspan=3)
	# ax2 = plt.subplot2grid((4,3),(2,0),colspan=3)
	ax3 = plt.subplot2grid((2,1),(1,0))

	ax0.set_xlabel(r'$x$, $c/\omega_{pl}$', fontsize=20)
	ax0.set_ylabel(r'$y$, $c/\omega_{pl}$', fontsize=20)
	# # ax1.set_xlabel(r'$x$, $c/\omega_{pl}$', fontsize=20)
	# ax1.set_ylabel(r'$y$, $c/\omega_{pl}$', fontsize=20)
	# ax2.set_xlabel(r'$x$, $c/\omega_{pl}$', fontsize=20)
	# ax2.set_ylabel(r'$y$, $c/\omega_{pl}$', fontsize=20)
	ax0.tick_params(axis='both', labelsize=20)
	# # ax0.tick_params(labelbottom='off')
	# ax1.tick_params(axis='both', labelsize=20)
	# # ax1.tick_params(labelbottom='off')
	# ax2.tick_params(axis='both', labelsize=20)
	# ax3.tick_params(axis='both', labelsize=20)
	ax0.set_aspect(1)
	# ax1.set_aspect(1)
	# ax2.set_aspect(1)

	n_colorbar = 5

	divider0 = make_axes_locatable(ax0)
	cax0 = divider0.append_axes("right", size="1%", pad=0.05)
	dens /= ppc0
	strm0 = ax0.pcolormesh(x, y, np.log(dens + 1.), cmap='plasma', vmin=np.log(1.01), vmax=np.log(maximum_dens + 1.01))
	cbar0 = plt.colorbar(strm0, cax = cax0)
	cbar0.set_ticks(np.linspace(np.log(1.01), np.log(maximum_dens + 1), n_colorbar))
	values = map(fmt, np.round(np.exp(np.linspace(np.log(0.01), np.log(maximum_dens), n_colorbar)),2))
	cbar0.set_ticklabels(values)
	cbar0.ax.yaxis.set_tick_params(pad=10)
	ax0.set_xlim(x.min(), x.max())
	ax0.set_ylim(-100 / (skin_depth / 10.), 100 / (skin_depth / 10.))
	cbar0.ax.tick_params(labelsize=20)

	# bx = getBField(root, step, 'bx', getSizes(root, step))
	# by = getBField(root, step, 'by', getSizes(root, step))
	# bz = getBField(root, step, 'bz', getSizes(root, step))
	# ax0.streamplot(x, y, by, bx, color = 'white', density = 0.6, linewidth = 1.2)

	# divider1 = make_axes_locatable(ax1)
	# cax1 = divider1.append_axes("right", size="1%", pad=0.05)
	# denbw /= ppc0
	# strm1 = ax1.pcolormesh(x, y, np.log(denbw + 1.), cmap='YlGnBu_r', vmin=np.log(1.01), vmax=np.log(maximum_denbw + 1.01))
	# cbar1 = plt.colorbar(strm1, cax = cax1)
	# cbar1.set_ticks(np.linspace(np.log(1.01), np.log(maximum_denbw + 1), n_colorbar))
	# values = map(fmt, np.round(np.exp(np.linspace(np.log(1.01), np.log(maximum_denbw + 1.01), n_colorbar)) - 1., 2))
	# cbar1.set_ticklabels(values)
	# cbar1.ax.yaxis.set_tick_params(pad=10)
	# ax1.set_xlim(x.min(), x.max())
	# ax1.set_ylim(-100 / (skin_depth / 10.), 100 / (skin_depth / 10.))
	# cbar1.ax.tick_params(labelsize=20)
	#
	# divider2 = make_axes_locatable(ax2)
	# cax2 = divider2.append_axes("right", size="1%", pad=0.05)
	# strm2 = ax2.pcolormesh(x, y, np.log(densph + 1.), cmap='afmhot', vmin=np.log(1.01), vmax=np.log(maximum_densph + 1.01))
	# cbar2 = plt.colorbar(strm2, cax = cax2)
	# cbar2.set_ticks(np.linspace(np.log(1.01), np.log(maximum_densph + 1), n_colorbar))
	# values = map(fmt, np.round(np.exp(np.linspace(np.log(1.01), np.log(maximum_densph + 1.01), n_colorbar)) - 1., 2))
	# cbar2.set_ticklabels(values)
	# cbar2.ax.yaxis.set_tick_params(pad=10)
	# ax2.set_xlim(x.min(), x.max())
	# ax2.set_ylim(-100 / (skin_depth / 10.), 100 / (skin_depth / 10.))
	# cbar2.ax.tick_params(labelsize=20)

	if len(photons_ch) > 0:
		cnts, bns = np.histogram(np.log(photons + 1), bins=200, weights = photons_ch)
		bns = average(bns)
		bns = np.exp(bns) - 1
		cnts = bns * cnts
		ax3.step(bns, cnts, color = 'black', label = 'photons')

	# if len(pairs) > 0:
	# 	cnts, bns = np.histogram(np.log(pairs), bins=200)
	# 	bns = average(bns)
	# 	bns = np.exp(bns)
	# 	cnts = bns * cnts
	# 	ax3.step(bns, cnts, color = 'red', label = r'pairs')

	# cnts, bns = np.histogram(np.log(particles), bins=200)
	# bns = average(bns)
	# bns = np.exp(bns)
	# cnts = cnts * bns
	# ax3.step(bns, cnts, color = 'blue', label = r'plasma')

	ax3.set_xscale('log')
	ax3.set_yscale('log', nonposy='clip')
	ax3.yaxis.tick_right()
	ax3.yaxis.set_label_position("right")

	ax3.legend(loc='upper center', ncol=3, fontsize=20)
	ax3.ticklabel_format(fontsize=20)

	ax3.set_xlim(0.1, 1000)
	ax3.set_ylim(0, 1e5)
	ax3.set_xlabel(r'$\gamma$', fontsize=20)
	ax3.set_ylabel(r'$\gamma~f(\gamma)$', fontsize=20)

	props = dict(boxstyle='square', facecolor='white', alpha=0.9, edgecolor='none')
	ax0.text(0.02, 0.95, r'photons', transform=ax0.transAxes, fontsize=25, verticalalignment='top', bbox=props)
	# ax1.text(0.02, 0.95, r'pairs', transform=ax1.transAxes, fontsize=25, verticalalignment='top', bbox=props)
	# ax2.text(0.02, 0.95, r'photons', transform=ax2.transAxes, fontsize=25, verticalalignment='top', bbox=props)

	hlp.printProgress(step - start, min(max_number, end) - start - 1, prefix = 'Progress:', barLength = 50)
	plt.tight_layout()
	plt.savefig(output_dir + "all_" + str(step+1).zfill(3) + ".png", dpi=100)

#plt.show()
