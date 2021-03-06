import glob 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm
import os
import matplotlib as mpl 
from scipy.stats import linregress
from matplotlib import rc_file
from sys import exit

rcFile = 'rc_files/paper.rc'  

nameToColorKeys = {
        "F145X":'#A93226',
        "M218X":'r',
        "F165X":'#E67E22',
        "F114X":'y',
        "N198X":'#2ECC71',
        "Y143X":'g',
        "N212X":'c',
        "D117X":'b',
        "Y92X":'m'
        }

molecList = [
        "Y92X",
        "F114X", 
        "D117X", 
        "Y143X", 
        "F145X", 
        "F165X", 
        "N198X",
        "N212X", 
        "M218X"]

figCols=3
figRows=3

rc_file(rcFile) 

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

datafiles = glob.glob('B_State/*/hbond_with_ca/geometry.xvg') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.10,hspace=0.35,left=0.12,right=0.98,top=0.93,bottom=0.1) 
fig.text(0.5,0.02, r"$d_{\rm{NH}}$ ($\rm{\AA}$)", ha='center', va='center') 
fig.text(0.01,0.5, r"$\theta_1$ (deg)", ha='left', va='center',rotation='vertical') 

xmin,xmax = 1.45, 2.45
ymin,ymax = 99,180 

for index,molec in enumerate(molecList) : 
    file="B_State/GFP_%s/hbond_with_ca/geometry.xvg"%molec
    data = np.genfromtxt(file,skip_header=23) 
    try :
        data[:,0] = data[:,0] / 1000 * 4
    except IndexError : 
        print "%s is empty. Skipping"%file      
        data = np.empty((0,4),int)  
        data = np.append(data,[[0,0,0,0]],axis=0) 

    for i in range(len(data[:,3]))  : 
        if data[i,3] > 180 : 
            data[i,3] -= 180 
    ax = axarr[index/figCols,index%figCols]

    xbins, ybins = np.arange(xmin, xmax, 0.01), np.arange(ymin,ymax,1)
    z, x, y = np.histogram2d(data[:,2],data[:,3],[xbins,ybins])

    im = ax.pcolor(x,y,z.T, cmap='plasma', vmin = 0, vmax = 10) 

    ## Overlay contour lines 
    x = data[:,2] ; y = data[:,3]
    counts,  xbins,  ybins  = np.histogram2d(x, y,range=((xmin,xmax),(ymin,ymax))) #,bins=(64,64)) 
    #extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

#    ax.imshow(heatmap) 
#    plt.close() 

    ax.contour(counts.transpose(),extent=[xmin,xmax,ymin,ymax],
       #ybins.min(),ybins.max()],linewidths=1,colors='black',
       linewidths=1,colors='black',
           linestyles='solid',levels=np.arange(-1,1500,125))

    ax.set_title(molec,color=nameToColorKeys[molec]) 
    ax.set_ylim([ymin,179])
    ax.set_xlim([xmin,xmax]) 

fig.subplots_adjust(right=0.88) 
cbar_ax = fig.add_axes([0.90, 0.15, 0.015, 0.7]) 
fig.text(0.98,0.5, r"Counts per bin", ha='center', va='center',rotation='vertical') 
fig.colorbar(im, cax=cbar_ax) 


fig.savefig('figures/Geometries_heat.png',format='png') 
plt.close() 

#
#datafiles = glob.glob('B_State/*/hbond_with_ca/geometry.xvg') 
#
#fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
#fig.subplots_adjust(wspace=0.15,hspace=0.25) 
#fig.text(0.5,0.02, r"$\theta_2$ (deg)",ha='center', va='center') 
#fig.text(0.05,0.5, r"$\theta_1$ (deg)", ha='center', va='center',rotation='vertical') 
#
#xmin,xmax = 120,180
#ymin,ymax = 99,180 
#
#for index,file in enumerate(datafiles) : 
#    data = np.genfromtxt(file,skip_header=23) 
#    try :
#        data[:,0] = data[:,0] / 1000 * 4
#    except IndexError : 
#        print "%s is empty. Skipping"%file      
#        data = np.empty((0,4),int)  
#        data = np.append(data,[[0,0,0,0]],axis=0) 
#
#    for i in range(len(data[:,3]))  : 
#        if data[i,3] > 180 : 
#            data[i,3] -= 180 
#    ax = axarr[index/figCols,index%figCols]
#
#    xbins, ybins = np.arange(xmin, xmax, 1.00), np.arange(ymin,ymax,1)
#    z, x, y = np.histogram2d(data[:,4],data[:,3],[xbins,ybins])
#
#    ax.pcolor(x,y,z.T, cmap='plasma', vmin = 0, vmax = 10) 
#
#    ## Overlay contour lines 
#    x = data[:,4] ; y = data[:,3]
#    counts,  xbins,  ybins  = np.histogram2d(x, y,range=((xmin,xmax),(ymin,ymax))) #,bins=(64,64)) 
#    #extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
#
##    ax.imshow(heatmap) 
##    plt.close() 
#
#    ax.contour(counts.transpose(),extent=[xmin,xmax,ymin,ymax],
#       #ybins.min(),ybins.max()],linewidths=1,colors='black',
#       linewidths=1,colors='black',
#           linestyles='solid',levels=np.arange(-1,1500,125))
#
#    ax.set_title(file.split('/')[1].split('_')[1]) 
#    ax.set_ylim([ymin,ymax])
#    ax.set_xlim([xmin,xmax]) 
#
#
#if len(datafiles)%figCols != 0 : 
#    for i in np.arange(figCols, len(datafiles)%figCols, -1) : 
#        fig.delaxes(axarr[figRows-1,i-1]) 
#
#fig.savefig('figures/Geometries_theta1_vs_theta2.png',format='png') 
#plt.close() 
#
#
datafiles = glob.glob('B_State/*/hbond_with_ca/nw_geometry.xvg') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.10,hspace=0.35,left=0.12,right=0.98,top=0.93,bottom=0.1) 
fig.text(0.5,0.02, r"$d_{\rm{NH}}$ ($\rm{\AA}$)", ha='center', va='center') 
fig.text(0.01,0.5, r"$\theta_1$ (deg)", ha='left', va='center',rotation='vertical') 

xmin,xmax = 1.45, 2.45
ymin,ymax = 99,180 

for index,molec in enumerate(molecList) : 
    file="B_State/GFP_%s/hbond_with_ca/nw_geometry.xvg"%molec
    data = np.genfromtxt(file,skip_header=23) 
    try :
        data[:,0] = data[:,0] / 1000 * 4
    except IndexError : 
        print "%s is empty. Skipping"%file      
        data = np.empty((0,4),int)  
        data = np.append(data,[[0,0,0,0]],axis=0) #
    for i in range(len(data[:,3]))  : 
        if data[i,3] > 180 : 
            data[i,3] -= 180 
    ax = axarr[index/figCols,index%figCols]

    xbins, ybins = np.arange(xmin, xmax, 0.01), np.arange(ymin,ymax,1)
    z, x, y = np.histogram2d(data[:,2],data[:,3],[xbins,ybins])

    im = ax.pcolor(x,y,z.T, cmap='plasma', vmin = 0, vmax = 10) 

    ## Overlay contour lines 
    x = data[:,2] ; y = data[:,3]
    counts,  xbins,  ybins  = np.histogram2d(x, y,range=((xmin,xmax),(ymin,ymax))) #,bins=(64,64)) 
    #extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

#    ax.imshow(heatmap) 
#    plt.close() 

    ax.contour(counts.transpose(),extent=[xmin,xmax,ymin,ymax],
       #ybins.min(),ybins.max()],linewidths=1,colors='black',
       linewidths=1,colors='black',
           linestyles='solid',levels=np.arange(-1,1500,125))

    ax.set_title(molec,color=nameToColorKeys[molec]) 
    ax.set_ylim([ymin,179])
    ax.set_xlim([xmin,xmax]) 

fig.subplots_adjust(right=0.88) 
cbar_ax = fig.add_axes([0.90, 0.15, 0.015, 0.7]) 
fig.text(0.98,0.5, r"Counts per bin", ha='center', va='center',rotation='vertical') 
fig.colorbar(im, cax=cbar_ax) 

fig.savefig('figures/Geometries_protein_heat.png',format='png') 
plt.close() 
#
#datafiles = glob.glob('B_State/*/hbond_with_ca/nw_geometry.xvg') 
#
#fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
#fig.subplots_adjust(wspace=0.15,hspace=0.25) 
#fig.text(0.5,0.02, r"$\theta_2$ (deg)",ha='center', va='center') 
#fig.text(0.05,0.5, r"$\theta_1$ (deg)", ha='center', va='center',rotation='vertical') 
#
#xmin,xmax = 120,180
#ymin,ymax = 99,180 
#
#for index,file in enumerate(datafiles) : 
#    data = np.genfromtxt(file,skip_header=23) 
#    try :
#        data[:,0] = data[:,0] / 1000 * 4
#    except IndexError : 
#        print "%s is empty. Skipping"%file      
#        data = np.empty((0,4),int)  
#        data = np.append(data,[[0,0,0,0]],axis=0) 
#
#    for i in range(len(data[:,3]))  : 
#        if data[i,3] > 180 : 
#            data[i,3] -= 180 
#    ax = axarr[index/figCols,index%figCols]
#
#    xbins, ybins = np.arange(xmin, xmax, 1.00), np.arange(ymin,ymax,1)
#    z, x, y = np.histogram2d(data[:,4],data[:,3],[xbins,ybins]) 
#    ax.pcolor(x,y,z.T, cmap='plasma', vmin = 0, vmax = 10) 
#
#    ## Overlay contour lines 
#    x = data[:,4] ; y = data[:,3]
#    counts,  xbins,  ybins  = np.histogram2d(x, y,range=((xmin,xmax),(ymin,ymax))) #,bins=(64,64)) 
#    #extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
#
##    ax.imshow(heatmap) 
##    plt.close() 
#
#    ax.contour(counts.transpose(),extent=[xmin,xmax,ymin,ymax],
#       #ybins.min(),ybins.max()],linewidths=1,colors='black',
#       linewidths=1,colors='black',
#           linestyles='solid',levels=np.arange(-1,1500,125))
#
#    ax.set_title(file.split('/')[1].split('_')[1]) 
#    ax.set_ylim([ymin,ymax])
#    ax.set_xlim([xmin,xmax]) 
#
#
#if len(datafiles)%figCols != 0 : 
#    for i in np.arange(figCols, len(datafiles)%figCols, -1) : 
#        fig.delaxes(axarr[figRows-1,i-1]) 
#
#fig.savefig('figures/Geometries_protein_theta1_vs_theta2.png',format='png') 
#plt.close() 
#
###
###
#####Repeat using theta2 instead of theta1
#datafiles = glob.glob('B_State/*/hbond_with_ca/geometry.xvg') 
#
#fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
#fig.subplots_adjust(wspace=0.15,hspace=0.25) 
#fig.text(0.5,0.02, r"$d_{\rm{NH}}$ ($\rm{\AA}$)", ha='center', va='center') 
#fig.text(0.05,0.5, r"$\theta_2$ (deg)", ha='center', va='center',rotation='vertical') 
#
#xmin,xmax = 1.45, 2.5
#ymin,ymax = 120,180 
#
#for index,file in enumerate(datafiles) : 
#    data = np.genfromtxt(file,skip_header=23) 
#    try :
#        data[:,0] = data[:,0] / 1000 * 4
#    except IndexError : 
#        print "%s is empty. Skipping"%file      
#        data = np.empty((0,5),int)  
#        data = np.append(data,[[0,0,0,0,0]],axis=0) 
#
#    for i in range(len(data[:,4]))  : 
#        if data[i,4] > 180 : 
#            data[i,4] -= 180 
#    ax = axarr[index/figCols,index%figCols]
#
#    xbins, ybins = np.arange(xmin, xmax, 0.01), np.arange(ymin,ymax,1)
#    z, x, y = np.histogram2d(data[:,2],data[:,4],[xbins,ybins])
#
#    ax.pcolor(x,y,z.T, cmap='plasma', vmin = 0, vmax = 10) 
#
#    ## Overlay contour lines 
#    x = data[:,2] ; y = data[:,4]
#    counts,  xbins,  ybins  = np.histogram2d(x, y,range=((xmin,xmax),(ymin,ymax))) #,bins=(64,64)) 
#    #extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
#
##    ax.imshow(heatmap) 
##    plt.close() 
#
#    ax.contour(counts.transpose(),extent=[xmin,xmax,ymin,ymax],
#       #ybins.min(),ybins.max()],linewidths=1,colors='black',
#       linewidths=1,colors='black',
#           linestyles='solid',levels=np.arange(-1,1500,125))
#
#    ax.set_title(file.split('/')[1].split('_')[1]) 
#    ax.set_ylim([ymin,ymax])
#    ax.set_xlim([xmin,xmax]) 
#
#
#if len(datafiles)%figCols != 0 : 
#    for i in np.arange(figCols, len(datafiles)%figCols, -1) : 
#        fig.delaxes(axarr[figRows-1,i-1]) 
#
#fig.savefig('figures/Geometries_theta2_heat.png',format='png') 
#plt.close() 
#
#
#datafiles = glob.glob('B_State/*/hbond_with_ca/nw_geometry.xvg') 
#
#fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
#fig.subplots_adjust(wspace=0.15,hspace=0.25) 
#fig.text(0.5,0.02, r"$d_{\rm{NH}}$ ($\AA$)", ha='center', va='center') 
#fig.text(0.05,0.5, r"$\theta_2$ (deg)", ha='center', va='center',rotation='vertical') 
#
#xmin,xmax = 1.45, 2.5
#ymin,ymax = 120,180 
#
#for index,file in enumerate(datafiles) : 
#    data = np.genfromtxt(file,skip_header=23) 
#    try :
#        data[:,0] = data[:,0] / 1000 * 4
#    except IndexError : 
#        print "%s is empty. Skipping"%file      
#        data = np.empty((0,5),int)  
#        data = np.append(data,[[0,0,0,0,0]],axis=0) 
#
#    for i in range(len(data[:,4]))  : 
#        if data[i,4] > 180 : 
#            data[i,4] -= 180 
#    ax = axarr[index/figCols,index%figCols]
#
#    xbins, ybins = np.arange(xmin, xmax, 0.01), np.arange(ymin,ymax,1)
#    z, x, y = np.histogram2d(data[:,2],data[:,4],[xbins,ybins])
#
#    ax.pcolor(x,y,z.T, cmap='plasma', vmin = 0, vmax = 10) 
#
#    ## Overlay contour lines 
#    x = data[:,2] ; y = data[:,4]
#    counts,  xbins,  ybins  = np.histogram2d(x, y,range=((xmin,xmax),(ymin,ymax))) #,bins=(64,64)) 
#    #extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
#
##    ax.imshow(heatmap) 
##    plt.close() 
#
#    ax.contour(counts.transpose(),extent=[xmin,xmax,ymin,ymax],
#       #ybins.min(),ybins.max()],linewidths=1,colors='black',
#       linewidths=1,colors='black',
#           linestyles='solid',levels=np.arange(-1,1500,125))
#
#    ax.set_title(file.split('/')[1].split('_')[1]) 
#    ax.set_ylim([ymin,ymax])
#    ax.set_xlim([xmin,xmax]) 
#
#
#if len(datafiles)%figCols != 0 : 
#    for i in np.arange(figCols, len(datafiles)%figCols, -1) : 
#        fig.delaxes(axarr[figRows-1,i-1]) 
#
#fig.savefig('figures/Geometries_theta2_protein_heat.png',format='png') 
#plt.close() 
#
#
###Plot histogrammed data
#
#datafiles = glob.glob('B_State/*/fit_hbond_with_ca/dist.poly') 
#
###Dist analysis
fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.10,hspace=0.35,left=0.20,right=0.98,top=0.93,bottom=0.1) 
fig.text(0.5,0.02, r"$d_{\rm{NH}}$ ($\rm{\AA}$)",ha='center', va='center') 
fig.text(0.04,0.5, r"Occurences", ha='center', va='center',rotation='vertical') 

figDens, axarrDens = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
figDens.subplots_adjust(wspace=0.10,hspace=0.35,left=0.18,right=0.98,top=0.93,bottom=0.1) 
figDens.text(0.5,0.04, r"$d_{\rm{NH}}$ ($\rm{\AA}$)", ha='center', va='center') 
figDens.text(0.04,0.5, r"Probability density (Occurences/$\rm{\AA}^3$)", ha='center', va='center',rotation='vertical') 

f2, ax2 = plt.subplots(1,1,figsize=(4.3,3) )
f2.subplots_adjust(left=0.20,bottom=0.13,right=0.95,top=0.95) 
f2.text(0.5,0.04, r"$d_{\rm{NH}}$ ($\rm{\AA}$)", ha='center', va='center') 
f2.text(0.06,0.5, r"$\tilde{\nu}$ (cm$^{-1}$)", ha='center', va='center',rotation='vertical') 

absMax = {} 
with open('Exp_data/sasa_abs_data.dat') as f : 
    for line in f : 
        if not line.startswith('#') : 
            key = line.split()[0] 
            value = float(line.split()[2]) 
            absMax[key] = value 


accumAvgTheta = [] 
accumAbsMax = [] 
for index, molec in enumerate(molecList) : 
    fileWhis = 'B_State/GFP_%s/fit_hbond_with_ca/dist.his'%(molec) 
    fileWpoly = 'B_State/GFP_%s/fit_hbond_with_ca/dist.poly'%(molec) 
    fileWgeom = 'B_State/GFP_%s/fit_hbond_with_ca/geometry.xvg'%(molec) 
    filePhis = 'B_State/GFP_%s/fit_hbond_with_ca/nw_dist.his'%(molec) 
    filePpoly = 'B_State/GFP_%s/fit_hbond_with_ca/nw_dist.poly'%(molec) 
    filePgeom = 'B_State/GFP_%s/fit_hbond_with_ca/nw_geometry.xvg'%(molec) 
    
    headlines = 0 
    with open(fileWgeom) as f :
        lines = f.readlines() 
        for line in lines : 
            if line.startswith('#') or line.startswith('@') : 
                 headlines += 1 
            else : 
                 break 

    dataWhis = np.genfromtxt(fileWhis) 
    dataWpoly = np.genfromtxt(fileWpoly) 
    dataWgeom = np.genfromtxt(fileWgeom,skip_header=headlines) 

    headlines = 0 
    with open(fileWgeom) as f :
        lines = f.readlines() 
        for line in lines : 
            if line.startswith('#') or line.startswith('@') : 
                 headlines += 1 
            else : 
                 break 

    dataPhis = np.genfromtxt(filePhis) 
    dataPpoly = np.genfromtxt(filePpoly) 
    dataPgeom = np.genfromtxt(filePgeom,skip_header=headlines) 

    dataWhis[:,1] *= len(dataWgeom[:,0])  ## The histograms are unfortunately normalized. To un-normalize, 
    dataWpoly[:,1] *= len(dataWgeom[:,0])  # we mutiply each bin times the magnitude of the histogram (here), 
                                           # and the binSize (below in for loop) 
    dataPhis[:,1] *= len(dataPgeom[:,0]) 
    dataPpoly[:,1] *= len(dataPgeom[:,0]) 

    distsW = dataWpoly[:,0] 
    probsW = dataWpoly[:,1]

    distsP = dataPpoly[:,0] 
    probsP = dataPpoly[:,1]

    volumesW = np.zeros_like(distsW) 
    for i in range(len(volumesW)) : 
        if not i == len(volumesW) - 1 : 
            binSize = distsW[i+1] - distsW[i]
        else : 
            binSize = distsW[i] - distsW[i-1]
        dataWhis[i,1] *= binSize
        dataWpoly[i,1] *= binSize
        volumesW[i] = 4/3 * np.pi * ((distsW[i]+binSize/2)**3 - (distsW[i]-binSize/2)**3)
    #print volumesW
    probsW = probsW / volumesW

    volumesP = np.zeros_like(distsP) 
    for i in range(len(volumesP)) : 
        if not i == len(volumesP) - 1 : 
            binSize = distsP[i+1] - distsP[i]
        else : 
            binSize = distsP[i] - distsP[i-1]
        dataPhis[i,1] *= binSize
        dataPpoly[i,1] *= binSize
        volumesP[i] = 4/3 * np.pi * ((distsP[i]+binSize/2)**3 - (distsP[i]-binSize/2)**3)
    print np.sum(dataPhis[:,1])
    #print volumesP
    probsP = probsP / volumesP

    avgDistW = 0
    avgDistP = 0
    for i in range(len(distsW)) : 
        avgDistW = avgDistW + probsW[i] * distsW[i]
    for i in range(len(distsP)) : 
        avgDistP = avgDistP + probsP[i] * distsP[i]

    avgDistTot = (avgDistW + avgDistP) / (sum(probsP)+sum(probsW))
    avgDistP = avgDistP / sum(probsP) 
    avgDistW = avgDistW / sum(probsW) 

    print molec, avgDistTot

    ax = axarr[index/figCols,index%figCols]
    axD = axarrDens[index/figCols,index%figCols]

    if True :
        try : 
            ax2.scatter(avgDistTot, absMax["GFP_"+molec],marker='P',label=molec,color=nameToColorKeys[molec],edgecolor='none',s=100)  
            accumAbsMax.append(absMax["GFP_"+molec]) 
            accumAvgTheta.append(avgDistTot) 

        except KeyError : 
            print "No key found for %s"%file.split('/')[1]
    else : 
        ax2.scatter(avgDistTot, absMax["GFP_"+molec],marker='P',label=molec,color=nameToColorKeys[molec],edgecolor='k',s=100,alpha=0.25)  

    ax2.axhline(2227.5,color='k',linestyle='--') 
    ax2.axhline(2235.9,color='#66CDFF',linestyle='--',zorder=1) 
    
    ax.scatter(dataWhis[:,0],dataWhis[:,1],color='b',linestyle='--',s=0.1)
    ax.scatter(dataPhis[:,0],dataPhis[:,1],color='g',linestyle='--',s=0.1)
    ax.plot(dataWpoly[:,0],dataWpoly[:,1],color='b',linestyle='--') 
    ax.plot(dataPpoly[:,0],dataPpoly[:,1],color='g',linestyle='--') 
    axD.plot(distsW,probsW,color='b',linestyle='-') 
    axD.plot(distsP,probsP,color='g',linestyle='-') 
    ax.axvline(avgDistW, color='b', linestyle='--') 
    ax.axvline(avgDistP, color='g', linestyle='--') 
    ax.axvline(avgDistTot, color='k', linestyle='--') 
    axD.axvline(avgDistW, color='b', linestyle='--') 
    axD.axvline(avgDistP, color='g', linestyle='--') 
    axD.axvline(avgDistTot, color='k', linestyle='--') 

    ax.set_title(molec,color=nameToColorKeys[molec]) 
    #ax.set_ylim(0,30000)
    ax.set_xlim(1.45,2.45) 

    axD.set_title(molec,color=nameToColorKeys[molec]) 
    #axD.set_ylim(0,100000) 
    axD.set_xlim(1.45,2.45) 

ax2.legend(loc=2,edgecolor='k',framealpha=1.0) 
ax2.set_xlim([1.97,2.26]) 

slope,intercept,r_value,p_value,std_error = linregress(accumAvgTheta, accumAbsMax)
print "r = %f, p = %f"%(r_value,p_value) 

xs = np.linspace(min(accumAvgTheta), max(accumAvgTheta) ) 
ys = slope * xs + intercept
ax2.plot(xs, ys, label="r = %0.3f"%r_value,color='k')

f2.text(0.80,0.89,r"r = %0.3f"%r_value) 

fig.savefig('figures/Geometries_dist.png',format='png') 
figDens.savefig('figures/Geometries_dist_density.png',format='png') 
f2.savefig('figures/abs_max_vs_max_dist.pdf',format='pdf') 
plt.close() 
#
#
###Theta 1 analyis
fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.10,hspace=0.35,left=0.15,right=0.98,top=0.93,bottom=0.1) 
fig.text(0.5,0.02, r"$\theta_1$ (deg)", ha='center', va='center') 
fig.text(0.04,0.5, r"Occurences", ha='center', va='center',rotation='vertical') 

figDens, axarrDens = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
figDens.subplots_adjust(wspace=0.10,hspace=0.35,left=0.15,right=0.98,top=0.93,bottom=0.1) 
figDens.text(0.5,0.02, r"$\theta_1$ (deg)", ha='center', va='center') 
figDens.text(0.03,0.5, r"Probability density (Occurences/$\rm{\AA}^3$)", ha='center', va='center',rotation='vertical') 

f2, ax2 = plt.subplots(1,1,figsize=(4.3,3) )
f2.subplots_adjust(left=0.20,bottom=0.13,right=0.95,top=0.95) 
f2.text(0.5,0.02, r"$\theta_1$ (deg)", ha='center', va='center') 
f2.text(0.06,0.5, r"$\tilde{\nu}$ (cm$^{-1}$)", ha='center', va='center',rotation='vertical') 

absMax = {} 
with open('Exp_data/sasa_abs_data.dat') as f : 
    for line in f : 
        if not line.startswith('#') : 
            key = line.split()[0] 
            value = float(line.split()[2]) 
            absMax[key] = value 


accumAvgTheta = [] 
accumAbsMax = [] 
for index, molec in enumerate(molecList) : 
    fileWhis = 'B_State/GFP_%s/fit_hbond_with_ca/theta1.his'%(molec) 
    fileWpoly = 'B_State/GFP_%s/fit_hbond_with_ca/theta1.poly'%(molec) 
    fileWgeom = 'B_State/GFP_%s/fit_hbond_with_ca/geometry.xvg'%(molec) 
    filePhis = 'B_State/GFP_%s/fit_hbond_with_ca/nw_theta1.his'%(molec) 
    filePpoly = 'B_State/GFP_%s/fit_hbond_with_ca/nw_theta1.poly'%(molec) 
    filePgeom = 'B_State/GFP_%s/fit_hbond_with_ca/nw_geometry.xvg'%(molec) 
    
    headlines = 0 
    with open(fileWgeom) as f :
        lines = f.readlines() 
        for line in lines : 
            if line.startswith('#') or line.startswith('@') : 
                 headlines += 1 
            else : 
                 break 

    dataWhis = np.genfromtxt(fileWhis) 
    dataWpoly = np.genfromtxt(fileWpoly) 
    dataWgeom = np.genfromtxt(fileWgeom,skip_header=headlines) 

    headlines = 0 
    with open(fileWgeom) as f :
        lines = f.readlines() 
        for line in lines : 
            if line.startswith('#') or line.startswith('@') : 
                 headlines += 1 
            else : 
                 break 

    dataPhis = np.genfromtxt(filePhis) 
    dataPpoly = np.genfromtxt(filePpoly) 
    dataPgeom = np.genfromtxt(filePgeom,skip_header=headlines) 

    dataWhis[:,1] *= len(dataWgeom[:,0])  ## The histograms are unfortunately normalized. To un-normalize, 
    dataWpoly[:,1] *= len(dataWgeom[:,0])  # we mutiply each bin times the magnitude of the histogram (here),  
                                           # and the binSize (below in for loop) 
    dataPhis[:,1] *= len(dataPgeom[:,0]) 
    dataPpoly[:,1] *= len(dataPgeom[:,0]) 

    anglesW = dataWpoly[:,0] 
    probsW = dataWpoly[:,1]

    anglesP = dataPpoly[:,0] 
    probsP = dataPpoly[:,1]

    avgAngle = np.average(anglesW,weights=probsW) 
    variance = np.average((anglesW- avgAngle)**2,weights=probsW) 
    print "%10s  %3.1f"%(molec,np.sqrt(variance) ) 

    r = 2.45 

    volumesW = np.zeros_like(anglesW) 
    for i in range(len(volumesW)) : 
        if not i == len(volumesW) - 1 : 
            binSize = anglesW[i+1] - anglesW[i]
        else : 
            binSize = anglesW[i] - anglesW[i-1]
        dataWhis[i,1] *= binSize
        dataWpoly[i,1] *= binSize
        volumesW[i] = 2*np.pi * r**3 / 3 * (-np.cos((anglesW[i]+binSize/2) * np.pi / 180)  + np.cos((anglesW[i]-binSize/2)* np.pi / 180.) )
    #print volumesW
    print np.sum(dataWhis[:,1]) 
    probsW = probsW / volumesW

    volumesP = np.zeros_like(anglesP) 
    for i in range(len(volumesP)) : 
        if not i == len(volumesP) - 1 : 
            binSize = anglesP[i+1] - anglesP[i]
        else : 
            binSize = anglesP[i] - anglesP[i-1]
        dataPhis[i,1] *= binSize
        dataPpoly[i,1] *= binSize
        volumesP[i] = 2*np.pi * r**3 / 3 * (-np.cos((anglesP[i]+binSize/2) * np.pi / 180)  + np.cos((anglesP[i]-binSize/2)* np.pi / 180.) )
    #print volumesP
    probsP = probsP / volumesP

    avgAngleW = 0
    avgAngleP = 0
    for i in range(len(anglesW)) : 
        avgAngleW = avgAngleW + probsW[i] * anglesW[i]
    for i in range(len(anglesP)) : 
        avgAngleP = avgAngleP + probsP[i] * anglesP[i]

    avgAngleTot = (avgAngleW + avgAngleP) / (sum(probsP)+sum(probsW))
    avgAngleP = avgAngleP / sum(probsP) 
    avgAngleW = avgAngleW / sum(probsW) 

    #print molec, avgAngleTot

    ax = axarr[index/figCols,index%figCols]
    axD = axarrDens[index/figCols,index%figCols]

    #if not molec == "Y92X" and not molec == "F145X" : 
    if True : 
        try : 
            ax2.scatter(avgAngleTot, absMax["GFP_"+molec],marker='P',label=molec,color=nameToColorKeys[molec],edgecolor='none',s=100,zorder=10)  
            accumAbsMax.append(absMax["GFP_"+molec]) 
            accumAvgTheta.append(avgAngleTot) 

        except KeyError : 
            print "No key found for %s"%file.split('/')[1]
    else : 
        ax2.scatter(avgAngleTot, absMax["GFP_"+molec],marker='P',label=molec,color=nameToColorKeys[molec],edgecolor='k',s=100,alpha=0.25)  

    ax2.axhline(2227.5,color='k',linestyle='--') 
    ax2.axhline(2235.9,color='#66CDFF',linestyle='--') 
    
    ax.scatter(dataWhis[:,0],dataWhis[:,1],color='b',edgecolor='none',s=2.5) 
    ax.scatter(dataPhis[:,0],dataPhis[:,1],color='g',edgecolor='none',s=2.5) 
    ax.plot(dataWpoly[:,0],dataWpoly[:,1],color='b',linestyle='--')#,linewidth=0.5) 
    ax.plot(dataPpoly[:,0],dataPpoly[:,1],color='g',linestyle='--')#,linewidth=0.5) 
    axD.plot(anglesW,probsW,color='b',linestyle='-') 
    axD.plot(anglesP,probsP,color='g',linestyle='-') 
#    ax.axvline(avgAngleW, color='b', linestyle='--') 
#    ax.axvline(avgAngleP, color='g', linestyle='--') 
#    ax.axvline(avgAngleTot, color='k', linestyle='--') 
    axD.axvline(avgAngleW, color='b', linestyle='--',zorder=1) 
    axD.axvline(avgAngleP, color='g', linestyle='--',zorder=1) 
    axD.axvline(avgAngleTot, color='k', linestyle='--',zorder=1) 

    ax.set_title(molec,color=nameToColorKeys[molec]) 
    ax.set_ylim(0,450  ) 
    ax.set_xlim(100,180) 

    axD.set_title(molec,color=nameToColorKeys[molec]) 
    axD.set_ylim(0,3500 ) 
    axD.set_xlim(100,180) 

ax2.legend(loc=4,edgecolor='k',framealpha=1.0) 

slope,intercept,r_value,p_value,std_error = linregress(accumAvgTheta, accumAbsMax)
print "r = %f, p = %f"%(r_value,p_value) 

xs = np.linspace(min(accumAvgTheta), max(accumAvgTheta) ) 
ys = slope * xs + intercept
ax2.plot(xs, ys, label="r = %0.3f"%r_value,color='k')

f2.text(0.23,0.90,r"r = %0.3f"%r_value) 

fig.savefig('figures/Geometries_angles.png',format='png') 
figDens.savefig('figures/Geometries_angles_density.png',format='png') 
f2.savefig('figures/abs_max_vs_max_theta.pdf',format='pdf') 
plt.close() 
#
#for field in ['total_field','rxn_field','coloumb_field'] :
#    fig, axarr = plt.subplots(3,3) 
#    fig.set_size_inches(18.5, 10.5)
#    fName = 'force_calc_APBS/'+field+'.txt'
#    with open(fName) as f : 
#        starkShifts = {} 
#        for line in f.readlines() : 
#            items = line.split() 
#            key = items[0] 
#            value = items[2]
#            starkShifts[key] = value 
#    
#    for i in range(0,9) : 
#        ax = axarr[i/figCols,i%figCols] 
#        gamma = i /12.0
#        for index, molec in enumerate(molecList) : 
#            x = accumAvgTheta[index] 
#            y = accumAbsMax[index] 
#            stark = float(starkShifts["GFP_"+molec] ) * gamma
#            y -= stark
#            ax.scatter(x,y,marker='D',label=molec,color=nameToColorKeys[molec]) 
#            ax.set_title(r"$\gamma$ = %2.2f"%gamma) 
#    fig.savefig('force_calc_APBS/angle+%s.png'%field,format='png')         
#
#datafiles = glob.glob('B_State/*/fit_hbond_with_ca/dist.poly') 
#
###Theta 2 analyis
#fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
#fig.subplots_adjust(wspace=0.10,hspace=0.35,left=0.15,right=0.98,top=0.93,bottom=0.1) 
#fig.text(0.5,0.02, r"$\theta_2$", ha='center', va='center') 
#fig.text(0.04,0.5, r"Occurences", ha='center', va='center',rotation='vertical') 
#
#figDens, axarrDens = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
#figDens.subplots_adjust(wspace=0.10,hspace=0.35,left=0.15,right=0.98,top=0.93,bottom=0.1) 
#figDens.text(0.5,0.04, r"$\theta_2$", ha='center', va='center') 
#figDens.text(0.04,0.5, r"$\rho_{\rm{dens}}$ Occurences/$\rm{\AA}^3$", ha='center', va='center',rotation='vertical') 
#
#f2, ax2 = plt.subplots(1,1) 
#f2.subplots_adjust(left=0.20,right=0.95,top=0.95) 
#f2.text(0.5,0.04, r"$\theta_2$", ha='center', va='center') 
#f2.text(0.06,0.5, r"$\tilde{\nu}$ (cm$^{-1}$)", ha='center', va='center',rotation='vertical') 
#
#absMax = {} 
#with open('Exp_data/sasa_abs_data.dat') as f : 
#    for line in f : 
#        if not line.startswith('#') : 
#            key = line.split()[0] 
#            value = float(line.split()[2]) 
#            absMax[key] = value 
#
#
#accumAvgTheta = [] 
#accumAbsMax = [] 
#for index, molec in enumerate(molecList) : 
#    fileWhis = 'B_State/GFP_%s/fit_hbond_with_ca/theta2.his'%(molec) 
#    fileWpoly = 'B_State/GFP_%s/fit_hbond_with_ca/theta2.poly'%(molec) 
#    fileWgeom = 'B_State/GFP_%s/fit_hbond_with_ca/geometry.xvg'%(molec) 
#    filePhis = 'B_State/GFP_%s/fit_hbond_with_ca/nw_theta2.his'%(molec) 
#    filePpoly = 'B_State/GFP_%s/fit_hbond_with_ca/nw_theta2.poly'%(molec) 
#    filePgeom = 'B_State/GFP_%s/fit_hbond_with_ca/nw_geometry.xvg'%(molec) 
#    
#    headlines = 0 
#    with open(fileWgeom) as f :
#        lines = f.readlines() 
#        for line in lines : 
#            if line.startswith('#') or line.startswith('@') : 
#                 headlines += 1 
#            else : 
#                 break 
#
#    dataWhis = np.genfromtxt(fileWhis) 
#    dataWpoly = np.genfromtxt(fileWpoly) 
#    dataWgeom = np.genfromtxt(fileWgeom,skip_header=headlines) 
#
#    headlines = 0 
#    with open(fileWgeom) as f :
#        lines = f.readlines() 
#        for line in lines : 
#            if line.startswith('#') or line.startswith('@') : 
#                 headlines += 1 
#            else : 
#                 break 
#
#    dataPhis = np.genfromtxt(filePhis) 
#    dataPpoly = np.genfromtxt(filePpoly) 
#    dataPgeom = np.genfromtxt(filePgeom,skip_header=headlines) 
#
#    dataWhis[:,1] *= len(dataWgeom[:,0]) 
#    dataWpoly[:,1] *= len(dataWgeom[:,0]) 
#
#    dataPhis[:,1] *= len(dataPgeom[:,0]) 
#    dataPpoly[:,1] *= len(dataPgeom[:,0]) 
#
#    anglesW = dataWpoly[:,0] 
#    probsW = dataWpoly[:,1]
#
#    anglesP = dataPpoly[:,0] 
#    probsP = dataPpoly[:,1]
#
#    r = 2.45 
##    volumes = 2*np.pi * r**3 / 3 * ( -np.cos((np.pi - (anglesW * np.pi / 180)) + 1 * np.pi / 180) + np.cos(((np.pi - anglesW* np.pi / 180.) )  ))
#    volumes = 1 ##volume element not dependent on theta2
#    probsW = probsW / volumes 
#
##    volumes = 2*np.pi * r**3 / 3 * ( -np.cos((np.pi - (anglesP * np.pi / 180)) + 1 * np.pi / 180) + np.cos(((np.pi - anglesP* np.pi / 180.) )  ))
#    volumes = 1 ##volume element not dependent on theta2
#    probsP = probsP / volumes 
#
#    avgAngleW = 0
#    avgAngleP = 0
#    for i in range(len(anglesW)) : 
#        avgAngleW = avgAngleW + probsW[i] * anglesW[i]
#    for i in range(len(anglesP)) : 
#        avgAngleP = avgAngleP + probsP[i] * anglesP[i]
#
#    avgAngleTot = (avgAngleW + avgAngleP) / (sum(probsP)+sum(probsW))
#    avgAngleP = avgAngleP / sum(probsP) 
#    avgAngleW = avgAngleW / sum(probsW) 
#
#    print molec, avgAngleTot
#
#    ax = axarr[index/figCols,index%figCols]
#    axD = axarrDens[index/figCols,index%figCols]
#
#    if not molec == "Y143X" and not molec == "Y92X" : 
#        try : 
#            ax2.scatter(avgAngleTot, absMax["GFP_"+molec],marker='P',label=molec,color=nameToColorKeys[molec],edgecolor='k',s=100)  
#            accumAbsMax.append(absMax["GFP_"+molec]) 
#            accumAvgTheta.append(avgAngleTot) 
#
#        except KeyError : 
#            print "No key found for %s"%file.split('/')[1]
#    else : 
#        ax2.scatter(avgAngleTot, absMax["GFP_"+molec],marker='P',label=molec,color=nameToColorKeys[molec],edgecolor='k',s=100,alpha=0.25)  
#
#    ax2.axhline(2227.5,color='k',linestyle='--') 
#    ax2.axhline(2235.9,color='b',linestyle='--') 
#    
#    ax.scatter(dataWhis[:,0],dataWhis[:,1],color='b',linestyle='--',s=0.1)
#    ax.scatter(dataPhis[:,0],dataPhis[:,1],color='g',linestyle='--',s=0.1)
#    ax.plot(dataWpoly[:,0],dataWpoly[:,1],color='b',linestyle='--') 
#    ax.plot(dataPpoly[:,0],dataPpoly[:,1],color='g',linestyle='--') 
#    axD.plot(anglesW,probsW,color='b',linestyle='-') 
#    axD.plot(anglesP,probsP,color='g',linestyle='-') 
#    ax.axvline(avgAngleW, color='b', linestyle='--') 
#    ax.axvline(avgAngleP, color='g', linestyle='--') 
#    ax.axvline(avgAngleTot, color='k', linestyle='--') 
#    axD.axvline(avgAngleW, color='b', linestyle='--') 
#    axD.axvline(avgAngleP, color='g', linestyle='--') 
#    axD.axvline(avgAngleTot, color='k', linestyle='--') 
#
#    ax.set_title(molec) 
#    ax.set_ylim(0,450  ) 
#    ax.set_xlim(100,180) 
#
#    axD.set_title(molec) 
#    axD.set_ylim(0,4900 ) 
#    axD.set_xlim(100,180) 
#
#slope,intercept,r_value,p_value,std_error = linregress(accumAvgTheta, accumAbsMax)
#print "r = %f, p = %f"%(r_value,p_value) 
#
#xs = np.linspace(min(accumAvgTheta), max(accumAvgTheta) ) 
#ys = slope * xs + intercept
#ax2.plot(xs, ys, label="r = %0.3f"%r_value)
#
#ax2.legend(loc=3) 
#fig.savefig('figures/Geometries_theta2.png',format='png') 
#figDens.savefig('figures/Geometries_theta2_density.png',format='png') 
#f2.savefig('figures/abs_max_vs_max_theta2.pdf',format='pdf') 
#plt.close() 
