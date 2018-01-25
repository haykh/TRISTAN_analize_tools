import h5py

from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.lines as lines
import matplotlib.patches as patches
import matplotlib.colors as colors
import matplotlib

from mpl_toolkits.axes_grid1 import make_axes_locatable

import matplotlib
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True

root1 = '/tigress/PERSEUS/hakobyan/tristan_test/output/'
root2 = '/tigress/PERSEUS/hakobyan/tristan_test2/output/'

param = "dens"
start = int(raw_input("Start: "))
end = int(raw_input("End: "))

ymin = 3
ymax = -3
xmin = 3000
xmax = 5000

rates1=[]
rates2=[]

for i in range(start, end + 1):
	print i
	f1=h5py.File(root1 + 'flds.tot.' + str(i+1).zfill(3),'r')
	#print(f1.keys())
	#exit()
	ez=f1['ez'].value

	ez=ez[0,:,:]
	ez=np.transpose(ez)
	ez=np.rot90(ez)
	ez = ez[ymin:ymax,int(len(ez[0])*0.5)]

	#ez = ez[ymin:ymax,1500]

	by=f1['by'].value

	by=by[0,:,:]
	by=np.transpose(by)
	by=np.rot90(by)
	bup = by[ymax,-3]
	rate = np.sum(ez) / (len(ez) * bup)
	rates1.append(rate)

	f1=h5py.File(root2 + 'flds.tot.' + str(i+1).zfill(3),'r')
	#print(f1.keys())
	#exit()
	ez=f1['ez'].value

	ez=ez[0,:,:]
	ez=np.transpose(ez)
	ez=np.rot90(ez)
	ez = ez[ymin:ymax,int(len(ez[0])*0.5)]

	#ez = ez[ymin:ymax,1500]

	by=f1['by'].value

	by=by[0,:,:]
	by=np.transpose(by)
	by=np.rot90(by)
	bup = by[ymax,-3]
	rate = np.sum(ez) / (len(ez) * bup)
	rates2.append(rate)

time = len(rates1) * 500. / (10. / 0.45)
times = np.linspace(0, time, len(rates1))

rates1 = np.array(rates1)
rates2 = np.array(rates2)
#rates += 1
fig, ax = plt.subplots()
ax.plot(times, rates1, color='red')
ax.plot(times, rates2, color='blue')
#ax.set_yscale("log", nonposx='clip')
ax.set_ylabel(r"$\gamma$")
ax.set_xlabel(r"$t,~\omega_p^{-1}$")
#ax.set_ylim(0, 0.4)
plt.savefig("reconnection/rate.png", dpi=100)
#plt.show()
	
