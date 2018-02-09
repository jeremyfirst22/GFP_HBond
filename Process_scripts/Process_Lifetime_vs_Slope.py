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
fig.text(0.5,0.04, "Time (ns)", ha='center', va='center') 
#fig.text(0.5,0.04, "Frames", ha='center', va='center') 

with open('fits/water_chosen_decay.xvg') as f: 
    lines = f.readlines() 
    fitKeys= {} 
    for line in lines : 
        splitline = line.split() 
        fitKeys[splitline[0]] = splitline[1:] 

with open('fits/protein_chosen_decay.xvg') as f: 
    lines = f.readlines() 
    protFitKeys= {} 
    for line in lines : 
        splitline = line.split() 
        protFitKeys[splitline[0]] = splitline[1:] 

xmax = 0 
for index,file in enumerate(datafiles) : 
    try : 
        data = np.genfromtxt(file,skip_header=24) 
        data[:,0] = data[:,0] / 1000 * 4
    except : 
        continue 

    if (np.max(data[:,0]) > xmax) : xmax = np.max(data[:,0]) 

    
    for i in np.arange(len(data)) : 
        if data[i,5] == 0 : 
            print "i = %i, data = %i"%(i,data[i,5]) 
            break 
    print i 
    x = np.geomspace(data[1,0], data[i,0], 100) 
    print data[1,0]
    print x[0] 

    ax = axarr[index/figCols,index%figCols]
    if len(fitKeys[file.split('/')[1]] ) == 4 : 
        print "Single exponential"
        a, A, r, adjR = fitKeys[file.split('/')[1]]
        a, A, r, adjR = float(a), float(A), float(r), float(adjR) 
        A = A * 1000 / 4
        y = a*np.exp(-1*A * x) 
    elif len(fitKeys[file.split('/')[1]] ) == 6 : 
        print "Double exponential"
        ax.text(0.005,1,"*",color='b',fontsize='large') 
        a, A, b, B, r, adjR = fitKeys[file.split('/')[1]]
        a, A, b, B, r, adjR = float(a), float(A), float(b), float(B), float(r), float(adjR) 
        A = A * 1000 / 4
        B = B * 1000 / 4
        y = a*np.exp(-1*A * x) + b* np.exp(-1* B* x) 

    ##ax.loglog(data[:,0],data[:,4],color='tab:orange',label='Nearby water') 
    ax.loglog(data[:,0],data[:,5],color='b',linestyle='-', label='HBonding water') 
    ax.loglog(x,y,color='k',linestyle='-.') 
#    ax.semilogy(data[:,0],data[:,5],color='b',linestyle='-', label='HBonding water') 
#    ax.semilogy(x,y,color='k',linestyle='-.') 
    ax.text(.2,  600, "r = %.3f"%(r) , fontsize=10,color='b') 

    for i in np.arange(len(data)) : 
        if data[i,6] == 0 : 
            break 
    x = np.geomspace(data[1,0], data[i,0], 100) 

    if len(protFitKeys[file.split('/')[1]] ) == 4 : 
        print "Single exponential"
        a, A, r, adjR = protFitKeys[file.split('/')[1]]
        a, A, r, adjR = float(a), float(A), float(r), float(adjR) 
        A = A * 1000 / 4
        y = a*np.exp(-1*A * x) 
    elif len(protFitKeys[file.split('/')[1]] ) == 6 : 
        print "Double exponential"
        ax.text(0.01,1,"*",color='g',fontsize='large') 
        a, A, b, B, r, adjR = protFitKeys[file.split('/')[1]]
        a, A, b, B, r, adjR = float(a), float(A), float(b), float(B), float(r), float(adjR) 
        A = A * 1000 / 4
        B = B * 1000 / 4
        y = a*np.exp(-1*A * x) + b* np.exp(-1* B* x) 

    print file.split('/')[1]
    if not (file.split('/')[1] == "GFP_N212X" or file.split('/')[1] == "GFP_F114X") : 
        ax.loglog(data[:,0],data[:,6],color='g',label='Hbonding protein') 
        ax.loglog(x,y,color='k',linestyle='-.') 
#        ax.semilogy(data[:,0],data[:,6],color='g',label='Hbonding protein') 
#        ax.semilogy(x,y,color='k',linestyle='-.') 
        ax.text(.2,  200, "r = %.3f"%(r) , fontsize=10,color='g') 

    ax.set_title(file.split('/')[1].split('_')[1]) 
    #ax.set_ylim([0.8,2000]) 

for index,file in enumerate(datafiles) : 
    ax = axarr[index/figCols,index%figCols]
    ax.set_xlim(0.004  ,5) 
    ax.set_ylim(0.4,2000) 
    
fig.savefig('figures/Lifetime_hbonds.pdf',format='pdf') 
plt.close() 

