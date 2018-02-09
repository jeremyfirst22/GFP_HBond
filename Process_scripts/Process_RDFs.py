import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 

figCols=3
figRows=3

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0,hspace=0.35) 
fig.text(0.5,0.04, r"Distance ($\rm{\AA}$)", ha='center', va='center') 
fig.text(0.05,0.5, r"P(r)", ha='center', va='center',rotation='vertical') 

datafiles = glob.glob('B_State/*/rdf/wat_cnf.xvg') 
index=0
for datafile in datafiles : 
    print index, index/figCols, index%figRows
    ax = axarr[index/figCols,index%figCols]

    print datafile

    headerlines = 0 
    with open(datafile) as f : 
        for line in f.readlines() : 
            if line.startswith('#') or line.startswith('@') : 
                headerlines += 1 
            else : break 

    try : 
        data = np.genfromtxt(datafile,skip_header=headerlines) 
    except IOError : 
        print "No file found for %s"%(solvent)  
        continue

    data[:,0] = data[:,0] * 10 # nm -> Angstroms 

    ax.plot(data[:,0],data[:,1]) 

    ax.set_title(datafile.split('/')[1].split('_')[1]) 
    ax.set_xlim(0,15) 
    #ax.set_ylim(10,50.0) 

    index +=1
fig.savefig('figures/rdfs.png',format='png') 
plt.close() 
