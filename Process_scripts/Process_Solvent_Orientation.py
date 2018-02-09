import glob 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm
import os

figCols=4
figRows=3

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

datafiles = glob.glob('B_State/*/sorient/sori.xvg') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0) 
fig.text(0.5,0.04, r"cos($\theta_1$)", ha='center', va='center') 
fig.text(0.04,0.50, r"Distribution", ha='center', va='center',rotation='vertical') 

for index,file in enumerate(datafiles) : 
    with open(file) as f: 
        filelines = f.readlines() 
        for line in filelines : 
            if line.startswith('@ subtitle ') : 
                words = line.split() 
                shellSize = words[5] 
                break 
    data = np.genfromtxt(file,skip_header=25) 

    ax = axarr[index/figCols,index%figCols]

    ax.plot(data[:,0],data[:,1]) 
    ax.text(0.8,3.9,shellSize,ha='right',va='top') 

    ax.set_title(file.split('/')[1]) 
    ax.set_xlim([-1.1,1.1]) 
    ax.set_ylim([0,4.0]) 

fig.savefig('figures/sori.pdf',format='pdf') 
plt.close() 

datafiles = glob.glob('B_State/*/sorient/snor.xvg') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0) 
fig.text(0.5,0.04, r"cos($\theta_2$)", ha='center', va='center') 
fig.text(0.04,0.50, r"Distribution", ha='center', va='center',rotation='vertical') 

for index,file in enumerate(datafiles) : 
    with open(file) as f: 
        filelines = f.readlines() 
        for line in filelines : 
            if line.startswith('@ subtitle ') : 
                words = line.split() 
                shellSize = words[5] 
                break 
    data = np.genfromtxt(file,skip_header=25) 

    ax = axarr[index/figCols,index%figCols]

    ax.plot(data[:,0],data[:,1]) 
    ax.text(0.8,3.9,shellSize,ha='right',va='top') 

    ax.set_title(file.split('/')[1]) 
    ax.set_xlim([0,1.1]) 
    ax.set_ylim([0,6.0]) 

fig.savefig('figures/snor.pdf',format='pdf') 
plt.close() 

datafiles = glob.glob('B_State/*/sorient/sord.xvg') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0) 
fig.text(0.5,0.04, r"r (nm)", ha='center', va='center') 
fig.text(0.04,0.50, r"Avg. cos($\theta_n$)", ha='center', va='center',rotation='vertical') 

for index,file in enumerate(datafiles) : 
    data = np.genfromtxt(file,skip_header=25) 
    #data[:,0] = data[:,0] / 1000 * 4

    ax = axarr[index/figCols,index%figCols]

    ax.plot(data[:,0],data[:,1]) 
    ax.plot(data[:,0],data[:,2]) 

    ax.set_title(file.split('/')[1]) 
    #ax.set_xlim([0,12]) 
    ax.set_ylim([-1, 1]) 

fig.savefig('figures/sord.pdf',format='pdf') 
plt.close() 

datafiles = glob.glob('B_State/*/sorient/scum.xvg') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0) 
fig.text(0.5,0.04, r"r (nm)", ha='center', va='center') 
fig.text(0.04,0.50, r"Cummulative cos($\theta_n$)", ha='center', va='center',rotation='vertical') 

for index,file in enumerate(datafiles) : 
    data = np.genfromtxt(file,skip_header=25) 
    #data[:,0] = data[:,0] / 1000 * 4

    ax = axarr[index/figCols,index%figCols]

    ax.plot(data[:,0],data[:,1]) 
    ax.plot(data[:,0],data[:,2]) 

    ax.set_title(file.split('/')[1]) 
    #ax.set_xlim([0,12]) 
    #ax.set_ylim([0,10]) 

fig.savefig('figures/scum.pdf',format='pdf') 
plt.close() 

datafiles = glob.glob('B_State/*/sorient/scount.xvg') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0) 
fig.text(0.5,0.04, "Time (ps)", ha='center', va='center') 

for index,file in enumerate(datafiles) : 
    data = np.genfromtxt(file,skip_header=25) 
    #data[:,0] = data[:,0] / 1000 * 4

    ax = axarr[index/figCols,index%figCols]

    ax.plot(data[:,0],data[:,1]) 

    ax.set_title(file.split('/')[1]) 
    #ax.set_xlim([0,12]) 
    #ax.set_ylim([0,10]) 

fig.savefig('figures/scount.pdf',format='pdf') 
plt.close() 
