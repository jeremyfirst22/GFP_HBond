import glob 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm
import os

figCols=4
figRows=3

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

for state in ['A','B'] : 
    datafiles = glob.glob('%s_State/*/HBond_nit/*.persistent.xvg'%state) 
    
    fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
    fig.subplots_adjust(wspace=0) 
    fig.text(0.5,0.04, "Time (ns)", ha='center', va='center') 

    for index,file in enumerate(datafiles) : 
        data = np.genfromtxt(file,skip_header=23) 
        data[:,0] = data[:,0] / 1000 * 4
    
        ax = axarr[index/figCols,index%figCols]
    
        ax.plot(data[:,0],data[:,4]) 

        ax.set_title(file.split('/')[1]) 
        ax.set_xlim([0,05]) 
        ax.set_ylim([0,10]) 
    
    fig.savefig('figures/Persistant_hbonds_%s.png'%state,format='png') 
    plt.close() 
    
    
