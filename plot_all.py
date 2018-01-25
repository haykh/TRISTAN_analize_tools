import h5py
from matplotlib.mlab import bivariate_normal
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors

from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

name = 'test_pairs_new'
title = 'synchrotron cooling\n' + r'$\sigma=50, \gamma_c=40, \gamma_{\rm rad}=200$'
#title = 'synchrotron cooling,\n' + r'$\sigma=50, \gamma_c=40, \gamma_{\rm rad}=200, \tau=0.01$'
root = '/tigress/PERSEUS/hakobyan/jul_aug_2017/test_runs/' + name + '/output/'
#root = 'data/'
#root = '/tigress/hakobyan/tristan_' + name + '/output/'
photQ = True
pairQ = True

def gammas (array):
	partUi = array["ui"].value
	partUe = array["ue"].value
	partVi = array["vi"].value
	partVe = array["ve"].value
	partWi = array["wi"].value
	partWe = array["we"].value

	newparti = np.sqrt(1 + partUi**2 + partVi**2 + partWi**2)
	newparte = np.sqrt(1 + partUe**2 + partVe**2 + partWe**2)
	return {
		'ion': newparti,
		'lec': newparte
	}

def getParticles(i, photonsQ, pairsQ, size):
	if (photonsQ):
		f1 = h5py.File(root + 'phot.tot.' + str(i+1).zfill(3),'r')
	f2 = h5py.File(root + 'prtl.tot.' + str(i+1).zfill(3),'r')
	#print(f2.keys())
	ions = gammas(f2)['ion']
	lecs = gammas(f2)['lec']
	xions = f2['xi'].value
    	xlecs = f2['xe'].value
	if (pairsQ):
		pari = ions[f2["indi"].value < 0]
		pare = lecs[f2["inde"].value < 0]
		#lecs = lecs[f2["inde"].value >= 0]
		#ions = ions[f2["indi"].value >= 0]
	# 4 is downsampling parameter
    	wid = 60 * 4
    	ions = ions[(xions >= size * 4 * 0.5 - wid) & (xions <= size * 4 * 0.5 + wid)]
   	lecs = lecs[(xlecs >= size * 4 * 0.5 - wid) & (xlecs <= size * 4 * 0.5 + wid)]
	particles = np.append(ions, lecs)
	if (photonsQ):
		photons=f1["chp"].value
		#up = np.array(f1["up"].value)
		#vp = np.array(f1["vp"].value)
		#wp = np.array(f1["wp"].value)
		#photons = np.sqrt(up**2+vp**2+wp**2)
	if (pairsQ):
		pairs = np.append(pari, pare)
	if (pairsQ):
		return {
			'phot': photons,
			'pair': pairs,
			'part': particles
		}
	elif (photonsQ):
		return {
			'phot': photons,
			'part': particles
		}
	else:
		return {
			'part': particles
		}

def getSpectra(i):
	f1 = h5py.File(root + 'spect.' + str(i+1).zfill(3), 'r')
	gamma = 1. + f1["gamma"].value
	spec = f1["spece"].value + f1["specp"].value
	# taking cut in the center
	num = len(spec[0]) / 2
	spec = spec[:,num]
	return {
		'spec': spec,
		'gamm': gamma
	}

class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

start = int(raw_input("Start: "))
end = int(raw_input("End: "))

def getField(i, param, size):
	xmin = int((size * 0.5) - 300)
	xmax = int((size * 0.5) + 300)
	#xmin = 3
	#xmax = -3
	ymin = 3
	ymax = -3
	f1=h5py.File(root + 'flds.tot.' + str(i+1).zfill(3),'r')
	#print(f1.keys())
	hh=f1[param].value
	hh1=hh[0,:,:]
	hh1=np.transpose(hh1)
	hh1=np.rot90(hh1)
	hh1 = hh1[ymin:ymax,xmin:xmax]
	if not pairQ:
		hh1 = np.rot90(hh1, 3)
	return hh1

def getBField(i, param,size):
	xmin = int((size * 0.5) - 300)
	xmax = int((size * 0.5) + 300)
	ymin = 3
	ymax = -3
	f2=h5py.File(root + 'flds.tot.' + str(i+1).zfill(3),'r')
	#print(f1.keys())
	hh=f2[param].value
	hh2=hh[0,:,:]
	if pairQ:
		hh2 = np.rot90(hh2, 2)
		hh2 = hh2[ymin:ymax,xmin:xmax]
	else:
		hh2=np.transpose(hh2)
		hh2 = np.rot90(hh2)
		hh2 = hh2[ymin:ymax,xmin:xmax]
		hh2 = np.rot90(hh2, 3)
	return hh2

def determineSize(i):
	f1=h5py.File(root + 'flds.tot.' + str(i+1).zfill(3),'r')
	hh=f1['densi'].value
	hh1=hh[0,:,:]
	hh1=np.transpose(hh1)
	hh1=np.rot90(hh1)
	hh1 = hh1[1:3,3:-3]
	return len(hh1[0])

def average(array):
	newarr = array[:-1]
	for i in range(len(array) - 1):
		newarr[i] = 0.5 * (array[i] + array[i + 1])
	return newarr

for step in range(start, end + 1):
	print step
	size = determineSize(step)

	dens = getField(step, 'dens', size)

	if (photQ):
		densph = getField(step, 'densph', size)
	if (pairQ):
		denbw = getField(step, 'denbw', size)
	#ez = getField(step, 'ez', size)

	if (photQ):
		photons = getParticles(step, photQ, pairQ, size)['phot']
	if (pairQ):
		pairs = getParticles(step, photQ, pairQ, size)['pair']
	particles = getParticles(step, photQ, pairQ, size)['part']

	# 4 is downsampling parameter, 10 is skin depth
	x = (np.arange(len(dens[0])) - max(np.arange(len(dens[0]))) * 0.5) * 4 / 10.
	y = (np.arange(len(dens)) - max(np.arange(len(dens))) * 0.5) * 4 / 10.
	x,y = np.meshgrid(x,y)

	if pairQ:
		fig = plt.figure(figsize=(15, 10))
	elif photQ:
		fig = plt.figure(figsize=(18, 15))
	else:
		fig = plt.figure(figsize=(18, 10))

	if pairQ:
		ax0 = plt.subplot2grid((3,3),(0,0),rowspan=2)
		ax1 = plt.subplot2grid((3,3),(0,1),rowspan=2)
		ax2 = plt.subplot2grid((3,3),(0,2),rowspan=2)
		ax4 = plt.subplot2grid((3,3),(2,0),colspan=3)

		ax1.set_title("pair plasma")
		ax2.set_title("photons")
	elif photQ:
		ax0 = plt.subplot2grid((3,2),(0,0),colspan=2)
		ax2 = plt.subplot2grid((3,2),(1,0),colspan=2)
		ax4 = plt.subplot2grid((3,2),(2,0),colspan=2)
		ax2.set_title("photons")
	else:
		ax0 = plt.subplot2grid((2,2),(0,0),colspan=2)
		ax4 = plt.subplot2grid((2,2),(1,0),colspan=2)

	# here 500 is the time between outputs, 10 - skin depth, 0.45 - light speed
	plt.suptitle(title + ', ' + r'$\omega_{\rm pl}^{-1}=$' + str(step * 500. / (10. / 0.45)), fontsize=20, y=0.99)
	ax0.set_title("plasma")
	ax4.set_title("spectrum")

	ax0.set_xlabel(r'$x$, $c/\omega_{pl}$', fontsize=15)
	ax0.set_ylabel(r'$y$, $c/\omega_{pl}$', fontsize=15)
	# ax0.set_aspect(1)
	# if photQ:
	# 	ax2.set_aspect(1)
	# if pairQ:
	# 	ax1.set_aspect(1)

	divider0 = make_axes_locatable(ax0)
	cax0 = divider0.append_axes("right", size="5%", pad=0.05)
	dens += 1
	strm0 = ax0.pcolormesh(x, y, np.log(dens))
	cbar0 = plt.colorbar(strm0, cax = cax0)
	cbar0.set_ticks(np.linspace(
			np.log(dens.min()),
			np.log(dens.max()),
			10)
		)
	cbar0.set_ticklabels(
		np.floor((np.exp(np.linspace(
			np.log(dens.min()),
			np.log(dens.max()),
			10))-1) / 5.)
	)
	ax0.set_xlim(x.min(), x.max())
	ax0.set_ylim(y.min(), y.max())

	bx = getBField(step, 'bx', size)
	by = getBField(step, 'by', size)
	bz = getBField(step, 'bz', size)

	if pairQ:
		ax0.streamplot(x, y, -bx, -by,
			color = 'white',
			density = 1.5, linewidth = 0.5)
	else:
		ax0.streamplot(x, y, by, bx,
			color = 'white',
			density = 1.5, linewidth = 0.5)

	if (photQ):
		densph += 1
		divider2 = make_axes_locatable(ax2)
		cax2 = divider2.append_axes("right", size="5%", pad=0.05)
		strm2 = ax2.pcolormesh(x, y, np.log(densph))
		cbar2 = plt.colorbar(strm2, cax = cax2)
		cbar2.set_ticks(np.linspace(
				np.log(densph.min()),
				np.log(densph.max()),
				10)
			)
		cbar2.set_ticklabels(
			np.floor(np.exp(np.linspace(
				np.log(densph.min()),
				np.log(densph.max()),
				10))-1)
		)
		ax2.set_xlim(x.min(), x.max())
		ax2.set_ylim(y.min(), y.max())
	if (pairQ):
		denbw += 1
		divider1 = make_axes_locatable(ax1)
		cax1 = divider1.append_axes("right", size="5%", pad=0.05)
		strm1 = ax1.pcolormesh(x, y, np.log(denbw))
		cbar1 = plt.colorbar(strm1, cax = cax1)
		cbar1.set_ticks(np.linspace(
				np.log(denbw.min()),
				np.log(denbw.max()),
				10)
			)
		cbar1.set_ticklabels(
			np.floor(np.exp(np.linspace(
				np.log(denbw.min()),
				np.log(denbw.max()),
				10))-1)
		)
		ax1.set_xlim(x.min(), x.max())
		ax1.set_ylim(y.min(), y.max())
		
	if (photQ):
		cnts, bns = np.histogram(np.log(photons + 1), bins=200)
		bns = average(bns)
		bns = np.exp(bns) - 1
		cnts = bns * cnts
		ax4.step(bns, cnts, color = 'black', label = 'photons')

	if (pairQ):
		cnts, bns = np.histogram(np.log(pairs), bins=200)
		bns = average(bns)
		bns = np.exp(bns)
		cnts = bns * cnts
		ax4.step(bns, cnts, color = 'red', label = 'pairs')
	cnts, bns = np.histogram(np.log(particles), bins=200)
	bns = average(bns)
	bns = np.exp(bns)
	cnts = cnts * bns
	ax4.step(bns, cnts, color = 'blue', label = 'plasma')
	ax4.set_xscale('log')
	ax4.set_yscale('log', nonposy='clip')
	ax4.yaxis.tick_right()
	ax4.yaxis.set_label_position("right")

	ax4.legend(loc='upper left', fontsize=20)
	ax4.ticklabel_format(fontsize=20)

	if (photQ):
		ax4.set_xlim(0.1, 1000)
	else:
		ax4.set_xlim(1., 1000)
	ax4.set_ylim(0, 1e5)
	ax4.set_xlabel(r'$\gamma$')
	ax4.set_ylabel(r'$\gamma~f(\gamma)$')
	plt.savefig("reconnection/test/" + name + '/' + "all_" + str(step+1).zfill(3) + ".png", dpi=100)

#plt.show()
