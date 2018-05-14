import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 
import matplotlib.cm as cm 

figCols=4
figRows=3

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0) 
fig.text(0.5,0.04, r"R$_g$ ($\AA$)", ha='center', va='center') 
fig.text(0.03,0.5, r"RMSD ($\AA$)", ha='center', va='center',rotation='vertical') 

try : 
    binSize = int(sys.argv[1] ) 
except IndexError : 
    print "Usage: %s < binSize >"%(sys.argv[0]) 
    print "Setting binSize to default of 100" 
    binSize = 100 

datafiles = glob.glob('B_State/*/rgyr_310/without_ter.xvg') 
index=0
for datafile in datafiles : 
    molec = datafile.split('/')[1]
    datafile2 = "B_State/%s/rmsd_310/without_ter.xvg"%molec
    print index, index/figCols, index%figRows
    ax = axarr[index/figCols,index%figCols]

    print datafile
    try : 
        skip_lines = 0 
        with open(datafile) as f : 
            for line in f.readlines() : 
                if line.startswith('#') or line.startswith('@') : 
                    skip_lines += 1
                else : 
                    break 
        data = np.genfromtxt(datafile,skip_header=skip_lines) 
    except IOError : 
        print "No gryate file found for %s"%(molec) 
        continue
    try : 
        skip_lines = 0 
        with open(datafile2) as f : 
            for line in f.readlines() : 
                if line.startswith('#') or line.startswith('@') : 
                    skip_lines += 1
                else : 
                    break 
        data2= np.genfromtxt(datafile2,skip_header=skip_lines) 
    except IOError : 
        print "No RMSD file found for %s"%(molec) 
        continue

    x = data[:,1] * 10 # nm -> Angstroms
    y = data2[:,1] * 10 
    print len(x), len(y) 
    assert len(x) == len(y) 

    while len(x) % binSize != 0 : 
        x = x[:-1]
        y = y[:-1]
    assert len(x) % binSize == 0 
    assert len(y) % binSize == 0 
    
    xs = np.mean(x.reshape(-1,binSize),axis=1) 
    ys = np.mean(y.reshape(-1,binSize),axis=1) 

    colors = cm.brg(np.linspace(0,1,len(ys)) ) 

    ax.plot(x,y,color='k',linewidth=0.1,alpha=0.3) 
    for x,y,c in zip(xs,ys,colors) : 
        ax.scatter(x,y,s=0.9,color=c,marker='o') 

    ax.set_title(molec.split('_')[1])  
    ax.set_xlim(16.5,17.25) 
    ax.set_ylim(0.5,2.5) 

    index +=1
fig.savefig('figures/rmsd_v_rgyr_310.png',format='png') 

