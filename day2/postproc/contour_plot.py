## run module load "python/3.7.11" before running this script

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18})

nx = 129; ny = 65
fname = '../q1/T_xy_129_065_0000.dat'
data = np.genfromtxt(fname)
x = data[:,0].reshape(129,65); xplot = x[:,0]
y = data[:,1].reshape(129,65); yplot = y[0,:]
Temp = data[:,2].reshape(129,65)
Temp_ex = data[:,2].reshape(129,65)

fig, ax = plt.subplots()
c = ax.contourf(xplot, yplot, Temp.T)
ax.set_aspect('equal')
plt.colorbar(c,fraction=0.025, pad=0.1)
ax.xaxis.set_label_text(r'$x$',fontsize=20)
ax.yaxis.set_label_text(r'$y$',fontsize=20)
plt.savefig('../q1/temp.png',format="png",dpi=300,bbox_inches='tight')

fig, ax = plt.subplots()
c = ax.contourf(xplot, yplot, Temp_ex.T)
ax.set_aspect('equal')
plt.colorbar(c,fraction=0.025, pad=0.1)
ax.xaxis.set_label_text(r'$x$',fontsize=20)
ax.yaxis.set_label_text(r'$y$',fontsize=20)
plt.savefig('../q1/temp_ex.png',format="png",dpi=300,bbox_inches='tight')
