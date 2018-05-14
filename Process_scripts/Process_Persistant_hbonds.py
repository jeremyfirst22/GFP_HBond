import glob 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm
import os
from sys import exit

figCols=3
figRows=3

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

datafiles = glob.glob('B_State/*/hbond_2/persistent.xvg') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0,hspace=0.30) 
fig.subplots_adjust(left=0.1,right=0.8,bottom=0.1,top=0.9) 
fig.text(0.5,0.04, "Time (ns)", ha='center', va='center') 

xmax = 0 
for index,file in enumerate(datafiles) : 
    try : 
        data = np.genfromtxt(file,skip_header=24) 
        data[:,0] = data[:,0] / 1000 * 4
    except : 
        continue 

    if (np.max(data[:,0]) > xmax) : xmax = np.max(data[:,0]) 

    ax = axarr[index/figCols,index%figCols]

#    ax.loglog(data[:,0],data[:,4],color='c',label='Nearby water') 
    ax.loglog(data[:,0],data[:,5],color='b', label='HBonding water') 
    ax.loglog(data[:,0],data[:,6],color='g',label='Hbonding protein') 

    ax.set_title(file.split('/')[1].split('_')[1]) 
    ax.set_ylim([0.8,2000]) 

for index,file in enumerate(datafiles) : 
    ax = axarr[index/figCols,index%figCols]
    ax.set_xlim(0.004,40 ) 
    

fig.legend(loc='center right', fontsize='x-small') 
fig.savefig('figures/Persistant_hbonds.pdf',format='pdf') 
plt.close() 

