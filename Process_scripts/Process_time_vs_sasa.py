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
fig.subplots_adjust(wspace=0) 
fig.text(0.5,0.04, r"Time (ns)", ha='center', va='center') 
fig.text(0.04,0.5, r"SASA ($\AA^2$)", ha='center', va='center',rotation='vertical') 

molList = glob.glob('B_State/*/sasa') 
index=0
for mol in molList : 
    print index, index/figCols, index%figRows
    ax = axarr[index/figCols,index%figCols]

    datafile = mol+'/cnf_area.xvg'
    print datafile
    try : 
        data1 = np.genfromtxt(datafile,skip_header=27) 
        data1[:,0] = data1[:,0] / 1000 
        ax.scatter(data1[:,0],data1[:,2],color='k',s=0.1) 
    except IOError : 
        print "No file found for %s"%(datafile)  
    except ValueError :
        print "Trying %s again without bottom line."%datafile
        data1 = np.genfromtxt(datafile,skip_header=27,skip_footer=1) 
        data1[:,0] = data1[:,0] / 1000 
        ax.scatter(data1[:,0],data1[:,2],color='k',s=0.1) 

    datafile = mol+'/nit_4_atoms_area.xvg'
    try : 
        data2 = np.genfromtxt(datafile,skip_header=27) 
        data2[:,0] = data2[:,0] / 1000 
        ax.scatter(data2[:,0],data2[:,2],color='g',s=0.1) 
    except IOError : 
        print "No file found for %s"%(datafile)  
    except ValueError :
        print "Trying %s again without bottom line."%datafile
        data2 = np.genfromtxt(datafile,skip_header=27,skip_footer=1) 
        data2[:,0] = data2[:,0] / 1000 
        ax.scatter(data2[:,0],data2[:,2],color='g',s=0.1) 

    datafile = mol+'/nh_ct_area.xvg'
    try : 
        data3 = np.genfromtxt(datafile,skip_header=27) 
        data3[:,0] = data3[:,0] / 1000 
        ax.scatter(data3[:,0],data3[:,2],color='b',s=0.1) 
    except IOError : 
        print "No file found for %s"%(datafile)  
    except ValueError :
        print "Trying %s again without bottom line."%datafile
        data3 = np.genfromtxt(datafile,skip_header=27,skip_footer=1) 
        data3[:,0] = data3[:,0] / 1000 
        ax.scatter(data3[:,0],data3[:,2],color='b',s=0.1) 


    ax.set_title(datafile.split('/')[1]) 
    ax.set_xlim(0,50) 
    ax.set_ylim(0,3)     

    index +=1

fig.savefig('figures/sasa_v_time.png',format='png') 
plt.close() 
