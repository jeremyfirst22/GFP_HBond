import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 

figCols=4
figRows=3

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0,hspace=0.35) 
fig.text(0.5,0.04, "Time (ns)", ha='center', va='center') 
fig.text(0.05,0.5, r"Distance to nearest image ($\AA$)", ha='center', va='center',rotation='vertical') 
datafiles = glob.glob('B_State/*/minimage/mindist.xvg') 
index=0
for datafile in datafiles : 
    print index, index/figCols, index%figRows
    ax = axarr[index/figCols,index%figCols]

    print datafile
    try : 
        data = np.genfromtxt(datafile,skip_header=27) 
    except IOError : 
        print "No file found for %s"%(solvent)  
        continue
    except ValueError :
        print "Trying %s again without bottom line."%datafile
        data = np.genfromtxt(datafile,skip_header=27,skip_footer=1) 


    data[:,0] = data[:,0] / 1000 
    data[:,1] = data[:,1] * 10 # nm -> Angstroms
    cutoff = np.zeros_like(data[:,0]) 
    for i in range(len(cutoff) ) : 
        cutoff[i] = 10.0 * 2 

    ax.scatter(data[:,0],data[:,1],color='b',s=0.1) 
    ax.plot(data[:,0],cutoff,'k--') 
    ax.set_title(datafile.split('/')[1].split('_')[1]) 
    ax.set_xlim(0,50) 
    ax.set_ylim(10,50.0) 

    index +=1
fig.savefig('figures/minimage.png',format='png') 
plt.close() 
