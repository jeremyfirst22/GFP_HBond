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

datafiles = glob.glob('B_State/*/hbond/hb_count.xvg') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0,hspace=0.3) 
fig.subplots_adjust(left=0.1,right=0.8, bottom=0.1,top=0.9) 
fig.text(0.5,0.04, "Time (ns)", ha='center', va='center') 
fig.text(0.04,0.5, "Counts (x1000)", rotation='vertical',ha='center', va='center') 

for index,file in enumerate(datafiles) : 
    data = np.genfromtxt(file,skip_header=25) 
    data[:,1] = data[:,1] /1000
    data[:,2] = data[:,2] /1000
    data[:,3] = data[:,3] /1000
    data[:,4] = data[:,5] /1000
    data[:,5] = data[:,5] /1000
    data[:,6] = data[:,6] /1000

    ax = axarr[index/figCols,index%figCols]
    
    width=0.25
    barwidth=0.23 
    #ax.bar(data[:,0], data[:,1],width ) 
    #for i in np.arange(0,6,1) : 
    #    ax.bar(data[:,0]+width*i, data[:,i+1],width ) 
    #for i in np.arange(0,6,1) : 
    #    ax.plot(data[:,0]+i, data[:,i+1] ) 
    ax.bar(data[:,0]-1*width,data[:,2],barwidth,color='b',label='Water') 
    ax.bar(data[:,0]+0*width,data[:,5],barwidth,color='g',label='Protein') 
    ax.bar(data[:,0]+1*width,data[:,6],barwidth,color='k',label='Protein or Water') 


    ax.set_title(file.split('/')[1]) 
    ax.set_xlim([0-2*width,3+2*width]) 
    ax.set_xticks(np.arange(0,4,1) ) 
    #ax.set_ylim([0,10000]) 
    
    
fig.legend(loc='center right',fontsize='x-small') 
fig.savefig('figures/Number_of_hbonds.pdf',format='pdf') 
plt.close() 

