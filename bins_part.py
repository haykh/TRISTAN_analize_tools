import h5py

from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.lines as lines
import matplotlib.patches as patches
import matplotlib.colors as colors

from mpl_toolkits.axes_grid1 import make_axes_locatable

root = '/tigress/hakobyan/tristan/output/'

for it in range(20):
	f1=h5py.File(root + 'prtl.tot.' + str(it+1).zfill(3),'r')

	#print f1.keys()
	#exit()

	xes = f1["xe"].value
	yes = f1["ye"].value

	print len(xes)

	size = 330

	dens = [[0 for x in range(size)] for y in range(size)]

	for i in range(len(xes)):
		dens[int(xes[i])][int(yes[i])] += 1

	dens = np.array(dens)
	dens = np.rot90(dens)
	
	Nbig = size
	Nsmall = size
	#dens = dens.reshape([Nsmall, Nbig//Nsmall, Nsmall, Nbig//Nsmall]).mean(3).mean(1)

	x = np.arange(Nsmall)
	y = np.arange(Nsmall)
	x,y = np.meshgrid(x,y)

	maxPair = max([max(sub) for sub in dens])
	minPair = min([min(sub) for sub in dens])

	#plt.imshow(hh1, vmin=0, vmax=3.5)
	fig, ax = plt.subplots()
	#divider = make_axes_locatable(ax)
	#cax = divider.append_axes("right", size="5%", pad=0.05)
	strm = ax.pcolormesh(x, y, dens)
	#cbar = plt.colorbar(strm, cax = cax)

	ax.set_xlim(x.min(), x.max())
	ax.set_ylim(y.min(), y.max())
	ax.set_aspect(1)

	plt.savefig("test/denp_" + str(it+1).zfill(3) + ".png", dpi=100)
