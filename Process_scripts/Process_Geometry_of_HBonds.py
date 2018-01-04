import glob 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm
import os
import matplotlib as mpl 

figCols=4
figRows=3

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

datafiles = glob.glob('B_State/*/hbond/geometry.xvg') 

#fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
#fig.subplots_adjust(wspace=0) 
#fig.text(0.5,0.04, "Time (ns)", ha='center', va='center') 
#fig.text(0.08,0.5, "HBond Angle", ha='center', va='center',rotation='vertical') 
#
#for index,file in enumerate(datafiles) : 
#    data = np.genfromtxt(file,skip_header=23) 
#    try :
#        data[:,0] = data[:,0] / 1000 * 4
#    except IndexError : 
#        print "%s is empty. Skipping"%file  
#        continue 
#
#    for i in range(len(data[:,2]))  : 
#        if data[i,3] > 180 : 
#            data[i,3] -= 180 
#    ax = axarr[index/figCols,index%figCols]
#
#    ax.scatter(data[:,0],data[:,3],s=0.1) 
#    ax.set_title(file.split('/')[1]) 
#    ax.set_xlim([-5,55])
#
#fig.savefig('figures/Geometries_angles.png',format='png') 
#plt.close() 
#
#
#fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
#fig.subplots_adjust(wspace=0) 
#fig.text(0.5,0.04, "Time (ns)", ha='center', va='center') 
#fig.text(0.08,0.5, r"HBond Length $\AA$", ha='center', va='center',rotation='vertical') 
#
#for index,file in enumerate(datafiles) : 
#    data = np.genfromtxt(file,skip_header=23) 
#    try :
#        data[:,0] = data[:,0] / 1000 * 4
#    except IndexError : 
#        print "%s is empty. Skipping"%file  
#        continue 
#
#    for i in range(len(data[:,3]))  : 
#        if data[i,3] > 180 : 
#            data[i,3] -= 180 
#    ax = axarr[index/figCols,index%figCols]
#
#    ax.scatter(data[:,0],data[:,2],s=0.1) 
#    ax.set_title(file.split('/')[1]) 
#    ax.set_xlim([-5,55])
#
#fig.savefig('figures/Geometries_length.png',format='png') 
#plt.close()
#
#fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
#fig.subplots_adjust(wspace=0) 
#fig.text(0.5,0.04, r"HBond Length ($\AA$)", ha='center', va='center') 
#fig.text(0.08,0.5, r"HBond Angle (deg)", ha='center', va='center',rotation='vertical') 
#
#for index,file in enumerate(datafiles) : 
#    data = np.genfromtxt(file,skip_header=23) 
#    try :
#        data[:,0] = data[:,0] / 1000 * 4
#    except IndexError : 
#        print "%s is empty. Skipping"%file  
#        continue 
#
#    for i in range(len(data[:,3]))  : 
#        if data[i,3] > 180 : 
#            data[i,3] -= 180 
#    ax = axarr[index/figCols,index%figCols]
#
#    print file.split('/')[1], "\t", len(data[:,3]) 
#    ax.scatter(data[:,2],data[:,3],s=0.05) 
#    ax.set_title(file.split('/')[1]) 
#    ax.set_ylim([90,190])
#    ax.set_xlim([1.4,2.5])
#    print file.split('/')[1],"\t",len(data[:,3]) 
#
#fig.savefig('figures/Geometries_2D.png',format='png') 
#plt.close() 
#
#fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
#fig.subplots_adjust(wspace=0) 
#fig.text(0.5,0.04, r"HBond Length ($\AA$)", ha='center', va='center') 
#fig.text(0.08,0.5, r"HBond Angle (deg)", ha='center', va='center',rotation='vertical') 
#
#for index,file in enumerate(datafiles) : 
#    data = np.genfromtxt(file,skip_header=23) 
#    try :
#        data[:,0] = data[:,0] / 1000 * 4
#    except IndexError : 
#        print "%s is empty. Skipping"%file  
#        continue 
#
#    for i in range(len(data[:,3]))  : 
#        if data[i,3] > 180 : 
#            data[i,3] -= 180 
#    ax = axarr[index/figCols,index%figCols]
#
#    x = data[:,2] ; y = data[:,3]
#    counts,  xbins,  ybins  = np.histogram2d(x, y) #,bins=(64,64)) 
#    #extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
#
##    ax.imshow(heatmap) 
##    plt.close() 
#
#    ax.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),
#       ybins.min(),ybins.max()],linewidths=3,colors='black',
#           linestyles='solid')
#
#    ax.set_title(file.split('/')[1]) 
#    ax.set_ylim([90,190])
#    ax.set_xlim([1.4,2.5])
#
#fig.savefig('figures/Geometries_contours.png',format='png') 
#plt.close() 

datafiles = glob.glob('B_State/*/hbond/geometry.xvg') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.15,hspace=0.25) 
fig.text(0.5,0.04, r"HBond Length ($\AA$)", ha='center', va='center') 
fig.text(0.08,0.5, r"HBond Angle (deg)", ha='center', va='center',rotation='vertical') 

xmin,xmax = 1.45, 2.5
ymin,ymax = 99,180 

for index,file in enumerate(datafiles) : 
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

    ax.pcolor(x,y,z.T, cmap='plasma', vmin = 0, vmax = 10) 

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

    ax.set_title(file.split('/')[1]) 
    ax.set_ylim([ymin,ymax])
    ax.set_xlim([xmin,xmax]) 


if len(datafiles)%figCols != 0 : 
    for i in np.arange(figCols, len(datafiles)%figCols, -1) : 
        fig.delaxes(axarr[figRows-1,i-1]) 

fig.savefig('figures/Geometries_heat.png',format='png') 
plt.close() 


datafiles = glob.glob('B_State/*/hbond/nw_geometry.xvg') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.15,hspace=0.25) 
fig.text(0.5,0.04, r"HBond Length ($\AA$)", ha='center', va='center') 
fig.text(0.08,0.5, r"HBond Angle (deg)", ha='center', va='center',rotation='vertical') 

xmin,xmax = 1.45, 2.5
ymin,ymax = 99,180 

for index,file in enumerate(datafiles) : 
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

    ax.pcolor(x,y,z.T, cmap='plasma', vmin = 0, vmax = 10) 

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

    ax.set_title(file.split('/')[1]) 
    ax.set_ylim([ymin,ymax])
    ax.set_xlim([xmin,xmax]) 


if len(datafiles)%figCols != 0 : 
    for i in np.arange(figCols, len(datafiles)%figCols, -1) : 
        fig.delaxes(axarr[figRows-1,i-1]) 

fig.savefig('figures/Geometries_protein_heat.png',format='png') 
plt.close() 


##Repeat using theta2 instead of theta1
datafiles = glob.glob('B_State/*/hbond/geometry.xvg') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.15,hspace=0.25) 
fig.text(0.5,0.04, r"HBond Length ($\AA$)", ha='center', va='center') 
fig.text(0.08,0.5, r"HBond Angle (deg)", ha='center', va='center',rotation='vertical') 

xmin,xmax = 1.45, 2.5
ymin,ymax = 120,180 

for index,file in enumerate(datafiles) : 
    data = np.genfromtxt(file,skip_header=23) 
    try :
        data[:,0] = data[:,0] / 1000 * 4
    except IndexError : 
        print "%s is empty. Skipping"%file      
        data = np.empty((0,5),int)  
        data = np.append(data,[[0,0,0,0,0]],axis=0) 

    for i in range(len(data[:,4]))  : 
        if data[i,4] > 180 : 
            data[i,4] -= 180 
    ax = axarr[index/figCols,index%figCols]

    xbins, ybins = np.arange(xmin, xmax, 0.01), np.arange(ymin,ymax,1)
    z, x, y = np.histogram2d(data[:,2],data[:,4],[xbins,ybins])

    ax.pcolor(x,y,z.T, cmap='plasma', vmin = 0, vmax = 10) 

    ## Overlay contour lines 
    x = data[:,2] ; y = data[:,4]
    counts,  xbins,  ybins  = np.histogram2d(x, y,range=((xmin,xmax),(ymin,ymax))) #,bins=(64,64)) 
    #extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

#    ax.imshow(heatmap) 
#    plt.close() 

    ax.contour(counts.transpose(),extent=[xmin,xmax,ymin,ymax],
       #ybins.min(),ybins.max()],linewidths=1,colors='black',
       linewidths=1,colors='black',
           linestyles='solid',levels=np.arange(-1,1500,125))

    ax.set_title(file.split('/')[1]) 
    ax.set_ylim([ymin,ymax])
    ax.set_xlim([xmin,xmax]) 


if len(datafiles)%figCols != 0 : 
    for i in np.arange(figCols, len(datafiles)%figCols, -1) : 
        fig.delaxes(axarr[figRows-1,i-1]) 

fig.savefig('figures/Geometries_theta2_heat.png',format='png') 
plt.close() 


datafiles = glob.glob('B_State/*/hbond/nw_geometry.xvg') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.15,hspace=0.25) 
fig.text(0.5,0.04, r"HBond Length ($\AA$)", ha='center', va='center') 
fig.text(0.08,0.5, r"HBond Angle (deg)", ha='center', va='center',rotation='vertical') 

xmin,xmax = 1.45, 2.5
ymin,ymax = 120,180 

for index,file in enumerate(datafiles) : 
    data = np.genfromtxt(file,skip_header=23) 
    try :
        data[:,0] = data[:,0] / 1000 * 4
    except IndexError : 
        print "%s is empty. Skipping"%file      
        data = np.empty((0,5),int)  
        data = np.append(data,[[0,0,0,0,0]],axis=0) 

    for i in range(len(data[:,4]))  : 
        if data[i,4] > 180 : 
            data[i,4] -= 180 
    ax = axarr[index/figCols,index%figCols]

    xbins, ybins = np.arange(xmin, xmax, 0.01), np.arange(ymin,ymax,1)
    z, x, y = np.histogram2d(data[:,2],data[:,4],[xbins,ybins])

    ax.pcolor(x,y,z.T, cmap='plasma', vmin = 0, vmax = 10) 

    ## Overlay contour lines 
    x = data[:,2] ; y = data[:,4]
    counts,  xbins,  ybins  = np.histogram2d(x, y,range=((xmin,xmax),(ymin,ymax))) #,bins=(64,64)) 
    #extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

#    ax.imshow(heatmap) 
#    plt.close() 

    ax.contour(counts.transpose(),extent=[xmin,xmax,ymin,ymax],
       #ybins.min(),ybins.max()],linewidths=1,colors='black',
       linewidths=1,colors='black',
           linestyles='solid',levels=np.arange(-1,1500,125))

    ax.set_title(file.split('/')[1]) 
    ax.set_ylim([ymin,ymax])
    ax.set_xlim([xmin,xmax]) 


if len(datafiles)%figCols != 0 : 
    for i in np.arange(figCols, len(datafiles)%figCols, -1) : 
        fig.delaxes(axarr[figRows-1,i-1]) 

fig.savefig('figures/Geometries_theta2_protein_heat.png',format='png') 
plt.close() 
