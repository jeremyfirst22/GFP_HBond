import glob 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm
import os
from sys import exit

figCols=4
figRows=3

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

datafiles = glob.glob('B_State/*/HBond/persistent.xvg') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0) 
fig.text(0.5,0.04, "Time (ns)", ha='center', va='center') 

xmax = 0 
for index,file in enumerate(datafiles) : 
    data = np.genfromtxt(file,skip_header=24) 
    data[:,0] = data[:,0] / 1000 * 4

    if (np.max(data[:,0]) > xmax) : xmax = np.max(data[:,0]) 

    ax = axarr[index/figCols,index%figCols]

#    ax.semilogy(data[:,0],data[:,4]) 
#    ax.semilogy(data[:,0],data[:,5]) 
#    ax.semilogy(data[:,0],data[:,6]) 

    ax.loglog(data[:,0],data[:,4]) 
    ax.loglog(data[:,0],data[:,5]) 
    ax.loglog(data[:,0],data[:,6]) 

    ax.set_title(file.split('/')[1]) 
    ax.set_xlim([0,2.5]) 
    ax.set_ylim([0.9,1000]) 

for index,file in enumerate(datafiles) : 
    ax = axarr[index/figCols,index%figCols]
    ax.set_xlim(0.001,1.0) 
    
fig.savefig('figures/Persistant_hbonds.pdf',format='pdf') 
plt.close() 

