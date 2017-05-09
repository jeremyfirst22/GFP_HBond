import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm

figCols=4
figRows=3

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

for state in ['A','B'] : 
    datafiles = glob.glob('%s_State/*/min_dist/cnf_water.xvg'%state) 
    
    fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
    fig.subplots_adjust(wspace=0) 
    fig.text(0.5,0.04, "Time (ns)", ha='center', va='center') 
    fig.text(0.08,0.5, r"Distance to nearest water ($\AA$)", ha='center', va='center',rotation='vertical') 
    
    for index,file in enumerate(datafiles) : 
        data = np.genfromtxt(file,skip_header=23) 
        data[:,0] = data[:,0] / 1000 * 4
        data[:,1] = data[:,1] * 10 # nm -> Angstroms
    
        ax = axarr[index/figCols,index%figCols]
    
        ax.scatter(data[:,0],data[:,1],s=0.1) 
        ax.set_title(file.split('/')[1]) 
        ax.set_xlim([-5,55])
        ax.set_ylim([1,8]) 
    
    fig.savefig('figures/Min_dist_%s.png'%state,format='png') 
    plt.close() 
    
    
