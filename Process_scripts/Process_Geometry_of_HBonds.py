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
    datafiles = glob.glob('%s_State/*/HBond_nit/*.geometry.xvg'%state) 
    
    fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
    fig.subplots_adjust(wspace=0) 
    fig.text(0.5,0.04, "Time (ns)", ha='center', va='center') 
    fig.text(0.08,0.5, "HBond Angle", ha='center', va='center',rotation='vertical') 
    
    for index,file in enumerate(datafiles) : 
        data = np.genfromtxt(file,skip_header=23) 
        data[:,0] = data[:,0] / 1000 * 4
    
        for i in range(len(data[:,2]))  : 
            if data[i,3] > 180 : 
                data[i,3] -= 180 
        ax = axarr[index/figCols,index%figCols]
    
        ax.scatter(data[:,0],data[:,3],s=0.1) 
        ax.set_title(file.split('/')[1]) 
        ax.set_xlim([-5,55])
    
    fig.savefig('figures/Geometries_angles_%s.png'%state,format='png') 
    plt.close() 
    
    
    fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
    fig.subplots_adjust(wspace=0) 
    fig.text(0.5,0.04, "Time (ns)", ha='center', va='center') 
    fig.text(0.08,0.5, r"HBond Length $\AA$", ha='center', va='center',rotation='vertical') 
    
    for index,file in enumerate(datafiles) : 
        data = np.genfromtxt(file,skip_header=23) 
        data[:,0] = data[:,0] / 1000 * 4
    
        for i in range(len(data[:,3]))  : 
            if data[i,3] > 180 : 
                data[i,3] -= 180 
        ax = axarr[index/figCols,index%figCols]
    
        ax.scatter(data[:,0],data[:,2],s=0.1) 
        ax.set_title(file.split('/')[1]) 
        ax.set_xlim([-5,55])
    
    fig.savefig('figures/Geometries_length_%s.png'%state,format='png') 
    plt.close()
    
    fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
    fig.subplots_adjust(wspace=0) 
    fig.text(0.5,0.04, r"HBond Length ($\AA$)", ha='center', va='center') 
    fig.text(0.08,0.5, r"HBond Angle (deg)", ha='center', va='center',rotation='vertical') 
    
    for index,file in enumerate(datafiles) : 
        data = np.genfromtxt(file,skip_header=23) 
        data[:,0] = data[:,0] / 1000 * 4
    
        for i in range(len(data[:,3]))  : 
            if data[i,3] > 180 : 
                data[i,3] -= 180 
        ax = axarr[index/figCols,index%figCols]
    
        print file.split('/')[1], "\t", len(data[:,3]) 
        ax.scatter(data[:,2],data[:,3],s=0.05) 
        ax.set_title(file.split('/')[1]) 
        ax.set_ylim([50,190])
        ax.set_xlim([1.4,2.5])
        print file.split('/')[1],"\t",len(data[:,3]) 
    
    fig.savefig('figures/Geometries_2D_%s.png'%state,format='png') 
    plt.close() 
    
    fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
    fig.subplots_adjust(wspace=0) 
    fig.text(0.5,0.04, r"HBond Length ($\AA$)", ha='center', va='center') 
    fig.text(0.08,0.5, r"HBond Angle (deg)", ha='center', va='center',rotation='vertical') 
    plt.close() 
    
    for index,file in enumerate(datafiles) : 
        data = np.genfromtxt(file,skip_header=23) 
        data[:,0] = data[:,0] / 1000 * 4
    
        for i in range(len(data[:,3]))  : 
            if data[i,3] > 180 : 
                data[i,3] -= 180 
        ax = axarr[index/figCols,index%figCols]
    
        x = data[:,2] ; y = data[:,3]
        counts,  xbins,  ybins  = np.histogram2d(x, y) #,bins=(64,64)) 
        #extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    
    #    ax.imshow(heatmap) 
    #    plt.close() 
    
        ax.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),
           ybins.min(),ybins.max()],linewidths=3,colors='black',
               linestyles='solid')
    
        ax.set_title(file.split('/')[1]) 
        ax.set_ylim([50,190])
        ax.set_xlim([1.4,2.5])
    
    fig.savefig('figures/Geometries_contour_%s.png'%state,format='png') 
