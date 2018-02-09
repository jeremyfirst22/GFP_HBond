import glob 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm
import os
import matplotlib as mpl 

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

figCols=3
figRows=3

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

datafiles = glob.glob('B_State/*/hbond/geometry.xvg') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.15,hspace=0.25) 
fig.text(0.5,0.04, r"$d_{\rm{NH}}$ ($\rm{\AA}$)", ha='center', va='center') 
fig.text(0.05,0.5, r"$\theta_1$ (deg)", ha='center', va='center',rotation='vertical') 

xmin,xmax = 1.45, 2.5
ymin,ymax = 99,180 

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
#    xbins, ybins = np.arange(xmin, xmax, 0.01), np.arange(ymin,ymax,1)
#    z, x, y = np.histogram2d(data[:,2],data[:,3],[xbins,ybins])
#
#    ax.pcolor(x,y,z.T, cmap='plasma', vmin = 0, vmax = 10) 
#
#    ## Overlay contour lines 
#    x = data[:,2] ; y = data[:,3]
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
#fig.savefig('figures/Geometries_heat.png',format='png') 
#plt.close() 
#
#datafiles = glob.glob('B_State/*/hbond/geometry.xvg') 
#
#fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
#fig.subplots_adjust(wspace=0.15,hspace=0.25) 
#fig.text(0.5,0.04, r"$\theta_2$ (deg)",ha='center', va='center') 
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
#datafiles = glob.glob('B_State/*/hbond/nw_geometry.xvg') 
#
#fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
#fig.subplots_adjust(wspace=0.15,hspace=0.25) 
#fig.text(0.5,0.04, r"$d_{\rm{NH}}$ ($\rm{\AA}$)", ha='center', va='center') 
#fig.text(0.05,0.5, r"$\theta_1$ (deg)", ha='center', va='center',rotation='vertical') 
#
#xmin,xmax = 1.45, 2.5
#ymin,ymax = 99,180 
#
#for index,file in enumerate(datafiles) : 
#    data = np.genfromtxt(file,skip_header=23) 
#    try :
#        data[:,0] = data[:,0] / 1000 * 4
#    except IndexError : 
#        print "%s is empty. Skipping"%file      
#        data = np.empty((0,4),int)  
#        data = np.append(data,[[0,0,0,0]],axis=0) #
#    for i in range(len(data[:,3]))  : 
#        if data[i,3] > 180 : 
#            data[i,3] -= 180 
#    ax = axarr[index/figCols,index%figCols]
#
#    xbins, ybins = np.arange(xmin, xmax, 0.01), np.arange(ymin,ymax,1)
#    z, x, y = np.histogram2d(data[:,2],data[:,3],[xbins,ybins])
#
#    ax.pcolor(x,y,z.T, cmap='plasma', vmin = 0, vmax = 10) 
#
#    ## Overlay contour lines 
#    x = data[:,2] ; y = data[:,3]
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
#fig.savefig('figures/Geometries_protein_heat.png',format='png') 
#plt.close() 
#
#datafiles = glob.glob('B_State/*/hbond/nw_geometry.xvg') 
#
#fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
#fig.subplots_adjust(wspace=0.15,hspace=0.25) 
#fig.text(0.5,0.04, r"$\theta_2$ (deg)",ha='center', va='center') 
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
#datafiles = glob.glob('B_State/*/hbond/geometry.xvg') 
#
#fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
#fig.subplots_adjust(wspace=0.15,hspace=0.25) 
#fig.text(0.5,0.04, r"$d_{\rm{NH}}$ ($\rm{\AA}$)", ha='center', va='center') 
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
#datafiles = glob.glob('B_State/*/hbond/nw_geometry.xvg') 
#
#fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
#fig.subplots_adjust(wspace=0.15,hspace=0.25) 
#fig.text(0.5,0.04, r"$d_{\rm{NH}}$ ($\AA$)", ha='center', va='center') 
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
datafiles = glob.glob('B_State/*/fit_hbond/dist.poly') 
molecList = [] 
for file in datafiles : 
    molecList.append(file.split('/')[1]) 
print molecList 
#
#fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
#fig.subplots_adjust(wspace=0) 
#fig.text(0.5,0.04, "HBond distance", ha='center', va='center') 
#fig.text(0.08,0.5, "N", ha='center', va='center',rotation='vertical') 
#
#for index, molec in enumerate(molecList) : 
#    file = 'B_State/%s/fit_hbond/dist.his'%(molec) 
#    file2 = 'B_State/%s/fit_hbond/dist.poly'%(molec) 
#    file3 = 'B_State/%s/fit_hbond/geometry.xvg'%(molec) 
#    
#    headlines = 0 
#    with open(file3) as f :
#        lines = f.readlines() 
#        for line in lines : 
#            if line.startswith('#') or line.startswith('@') : 
#                 headlines += 1 
#            else : 
#                 break 
#
#    data = np.genfromtxt(file) 
#    data2 = np.genfromtxt(file2) 
#    data3 = np.genfromtxt(file3,skip_header=headlines) 
#
#    data[:,1] *= len(data3[:,0]) 
#    data2[:,1] *= len(data3[:,0]) 
#
#    ax = axarr[index/figCols,index%figCols]
#
#    ax.plot(data[:,0],data[:,1])
#    ax.plot(data2[:,0],data2[:,1])
#    ax.set_title(file.split('/')[1]) 
#    ax.set_ylim(0,40000) 
#    ax.set_xlim(1.5,2.5) 
#
#fig.savefig('figures/Geometries_dist.png',format='png') 
#plt.close() 
#
#fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
#fig.subplots_adjust(wspace=0) 
#fig.text(0.5,0.04, "HBond distance", ha='center', va='center') 
#fig.text(0.08,0.5, "N", ha='center', va='center',rotation='vertical') 
#
#for index, molec in enumerate(molecList) : 
#    file = 'B_State/%s/fit_hbond/nw_dist.his'%(molec) 
#    file2 = 'B_State/%s/fit_hbond/nw_dist.poly'%(molec) 
#    file3 = 'B_State/%s/fit_hbond/nw_geometry.xvg'%(molec) 
#    
#    headlines = 0 
#    with open(file3) as f :
#        lines = f.readlines() 
#        for line in lines : 
#            if line.startswith('#') or line.startswith('@') : 
#                 headlines += 1 
#            else : 
#                 break 
#
#    data = np.genfromtxt(file) 
#    data2 = np.genfromtxt(file2) 
#    data3 = np.genfromtxt(file3,skip_header=headlines) 
#
#    data[:,1] *= len(data3[:,0]) 
#    data2[:,1] *= len(data3[:,0]) 
#
#    ax = axarr[index/figCols,index%figCols]
#
#    ax.plot(data[:,0],data[:,1])
#    ax.plot(data2[:,0],data2[:,1])
#    ax.set_title(file.split('/')[1]) 
#    ax.set_ylim(0,40000) 
#    ax.set_xlim(1.5,2.5) 
#
#fig.savefig('figures/Geometries_nw_dist.png',format='png') 
#plt.close() 
#
#
fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0,hspace=0.25) 
fig.text(0.5,0.04, r"$\theta_1$", ha='center', va='center') 
fig.text(0.04,0.5, r"Occurances", ha='center', va='center',rotation='vertical') 

figDens, axarrDens = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
figDens.subplots_adjust(wspace=0,hspace=0.25) 
figDens.text(0.5,0.04, r"$\theta_1$", ha='center', va='center') 
figDens.text(0.04,0.5, r"$\rho_{\rm{dens}}$ Occurances/$\rm{\AA}^3$", ha='center', va='center',rotation='vertical') 

f2, ax2 = plt.subplots(1,1) 
f2.subplots_adjust(left=0.15) 
f2.text(0.5,0.04, r"$\theta_1$", ha='center', va='center') 
f2.text(0.04,0.5, r"$\tilde{\nu}$ (cm$^{-1}$", ha='center', va='center',rotation='vertical') 

absMax = {} 
with open('Exp_data/sasa_abs_data.dat') as f : 
    for line in f : 
        if not line.startswith('#') : 
            key = line.split()[0] 
            value = float(line.split()[2]) 
            absMax[key] = value 


for index, molec in enumerate(molecList) : 
    fileWhis = 'B_State/%s/fit_hbond/theta1.his'%(molec) 
    fileWpoly = 'B_State/%s/fit_hbond/theta1.poly'%(molec) 
    fileWgeom = 'B_State/%s/fit_hbond/geometry.xvg'%(molec) 
    filePhis = 'B_State/%s/fit_hbond/nw_theta1.his'%(molec) 
    filePpoly = 'B_State/%s/fit_hbond/nw_theta1.poly'%(molec) 
    filePgeom = 'B_State/%s/fit_hbond/nw_geometry.xvg'%(molec) 
    
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

    dataWhis[:,1] *= len(dataWgeom[:,0]) 
    dataWpoly[:,1] *= len(dataWgeom[:,0]) 

    dataPhis[:,1] *= len(dataPgeom[:,0]) 
    dataPpoly[:,1] *= len(dataPgeom[:,0]) 

    anglesW = dataWpoly[:,0] 
    probsW = dataWpoly[:,1]

    anglesP = dataPpoly[:,0] 
    probsP = dataPpoly[:,1]

    r = 2.45 
    volumes = np.pi * r**3 / 3 * ( -np.cos((np.pi - (anglesW * np.pi / 180)) + 1 * np.pi / 180) + np.cos(((np.pi - anglesW* np.pi / 180.) )  ))
    probsW = probsW / volumes 

    volumes = np.pi * r**3 / 3 * ( -np.cos((np.pi - (anglesP * np.pi / 180)) + 1 * np.pi / 180) + np.cos(((np.pi - anglesP* np.pi / 180.) )  ))
    probsP = probsP / volumes 

    avgAngleW = 0
    avgAngleP = 0
    for i in range(len(anglesW)) : 
        avgAngleW = avgAngleW + probsW[i] * anglesW[i]
    for i in range(len(anglesP)) : 
        avgAngleP = avgAngleP + probsP[i] * anglesP[i]

    avgAngleTot = (avgAngleW + avgAngleP) / (sum(probsP)+sum(probsW))
    avgAngleP = avgAngleP / sum(probsP) 
    avgAngleW = avgAngleW / sum(probsW) 

    print molec, avgAngleTot

    ax = axarr[index/figCols,index%figCols]
    axD = axarrDens[index/figCols,index%figCols]

    name = molec.split('_')[1]
    if not name == "Y143X" : 
        try : 
            ax2.scatter(avgAngleTot, absMax["GFP_"+name],marker='P',label=name,color=nameToColorKeys[name],edgecolor='k',s=100)  
            #ax2.scatter(avgAngleW, absMax["GFP_"+name],marker='o',color=nameToColorKeys[name],alpha=0.5)   
            #ax2.scatter(avgAngleP, absMax["GFP_"+name],marker='D',color=nameToColorKeys[name],alpha=0.5)   
        except KeyError : 
            print "No key found for %s"%file.split('/')[1]

    #ax.plot(data[:,0],data[:,1], color='b')
    ax.scatter(dataWhis[:,0],dataWhis[:,1],color='b',linestyle='--',s=0.1)
    ax.scatter(dataPhis[:,0],dataPhis[:,1],color='g',linestyle='--',s=0.1)
    ax.plot(dataWpoly[:,0],dataWpoly[:,1],color='b',linestyle='--') 
    ax.plot(dataPpoly[:,0],dataPpoly[:,1],color='g',linestyle='--') 
    axD.plot(anglesW,probsW,color='b',linestyle='-') 
    axD.plot(anglesP,probsP,color='g',linestyle='-') 
    ax.axvline(avgAngleW, color='b', linestyle='--') 
    ax.axvline(avgAngleP, color='g', linestyle='--') 
    ax.axvline(avgAngleTot, color='k', linestyle='--') 
    axD.axvline(avgAngleW, color='b', linestyle='--') 
    axD.axvline(avgAngleP, color='g', linestyle='--') 
    axD.axvline(avgAngleTot, color='k', linestyle='--') 

    ax.set_title(name) 
    ax.set_ylim(0,450  ) 
    ax.set_xlim(100,180) 

    axD.set_title(name) 
    axD.set_ylim(0,4900 ) 
    axD.set_xlim(100,180) 

ax2.legend(loc=4) 
fig.savefig('figures/Geometries_angles.png',format='png') 
figDens.savefig('figures/Geometries_angles_density.png',format='png') 
f2.savefig('figures/abs_max_vs_max_theta.png',format='png') 
plt.close() 

#fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
#fig.subplots_adjust(wspace=0) 
#fig.text(0.5,0.04, "Theta1", ha='center', va='center') 
#fig.text(0.08,0.5, "N", ha='center', va='center',rotation='vertical') 
#
##f2, ax2 = plt.subplots(1,1) 
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
#for index, molec in enumerate(molecList) : 
#    file = 'B_State/%s/fit_hbond/nw_theta1.his'%(molec) 
#    file2 = 'B_State/%s/fit_hbond/nw_theta1.poly'%(molec) 
#    file3 = 'B_State/%s/fit_hbond/nw_geometry.xvg'%(molec) 
#    
#    headlines = 0 
#    with open(file3) as f :
#        lines = f.readlines() 
#        for line in lines : 
#            if line.startswith('#') or line.startswith('@') : 
#                 headlines += 1 
#            else : 
#                 break 
#
#    data = np.genfromtxt(file) 
#    data2 = np.genfromtxt(file2) 
#    data3 = np.genfromtxt(file3,skip_header=headlines) 
#
#    data[:,1] *= len(data3[:,0]) 
#    data2[:,1] *= len(data3[:,0]) 
#
#    ax = axarr[index/figCols,index%figCols]
#
#    angles = data[:,0] 
#    probs = data2[:,1]
#
#    r = 2.45 
#    volumes = np.pi * r**3 / 3 * ( -np.cos((np.pi - (angles * np.pi / 180)) + 1 * np.pi / 180) + np.cos(((np.pi - angles* np.pi / 180.) )  )) 
#
#    probs = probs / volumes 
#
#    avgAngle = 0
#    for i in range(len(angles)) : 
#        avgAngle = avgAngle + probs[i] * angles[i]
#    avgAngle = avgAngle / sum(probs) 
#
#    imax, pmax = 0, 0 
#    for i in range(len(probs) ) : 
#        if probs[i] > pmax : 
#            imax, pmax = i, probs[i] 
#    print file.split('/')[1],angles[imax], pmax, avgAngle
#
#    #if not (file.split('/')[1] == 'GFP_Y92X' or file.split('/')[1] == 'GFP_Y143X' ) : 
#    molec = file.split('/')[1].split('_')[1]
#    if molec == "Y92X" : #or molec == "Y143X" : 
##    if True : 
#        try : 
#            ax2.scatter(avgAngle, absMax[file.split('/')[1]],label=file.split('/')[1].split('_')[1],marker='D')   
#        except KeyError : 
#            print "No key found for %s"%file.split('/')[1]
#    #elif molec == "F145X" : 
#    #    ax2.scatter(np.average([106.508741892,130.321072301]), absMax[file.split('/')[1]],label=molec,marker='+') 
#    elif molec == "F165X" : 
#        ax2.scatter(np.average([147.716857708,123.60041311]), absMax[file.split('/')[1]],label=molec,marker='+') 
#
#    #ax.plot(data[:,0],data[:,1],color='g')
#    ax.plot(data2[:,0],data2[:,1],color='g')
#    ax.axvline(avgAngle, color='k', linestyle='--') 
#    ax.plot(angles, probs,color='#2ECC71',linestyle='--') 
#
#    ax.set_title(file.split('/')[1]) 
#    #ax.set_ylim(0,450  ) 
#    #ax.set_xlim(100,180) 
#
#ax2.legend() 
#fig.savefig('figures/Geometries_nw_angles.png',format='png') 
#f2.savefig('figures/nw_abs_max_vs_max_theta.png',format='png') 
#plt.close() 

##
##
##
#fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
#fig.subplots_adjust(wspace=0) 
#fig.text(0.5,0.04, "Theta2", ha='center', va='center') 
#fig.text(0.08,0.5, "N", ha='center', va='center',rotation='vertical') 
#
#for index, molec in enumerate(molecList) : 
#    file = 'B_State/%s/fit_hbond/theta2.his'%(molec) 
#    file2 = 'B_State/%s/fit_hbond/theta2.poly'%(molec) 
#    file3 = 'B_State/%s/fit_hbond/geometry.xvg'%(molec) 
#    
#    headlines = 0 
#    with open(file3) as f :
#        lines = f.readlines() 
#        for line in lines : 
#            if line.startswith('#') or line.startswith('@') : 
#                 headlines += 1 
#            else : 
#                 break 
#
#    data = np.genfromtxt(file) 
#    data2 = np.genfromtxt(file2) 
#    data3 = np.genfromtxt(file3,skip_header=headlines) 
#
#    data[:,1] *= len(data3[:,0]) 
#    data2[:,1] *= len(data3[:,0]) 
#
#    ax = axarr[index/figCols,index%figCols]
#
#    ax.plot(data[:,0],data[:,1])
#    ax.plot(data2[:,0],data2[:,1])
#    ax.set_title(file.split('/')[1]) 
#    ax.set_ylim(0,450  ) 
#    ax.set_xlim(100,180) 
#
#fig.savefig('figures/Geometries_theta2_angles.png',format='png') 
#plt.close() 
##
#fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
#fig.subplots_adjust(wspace=0) 
#fig.text(0.5,0.04, "Theta2", ha='center', va='center') 
#fig.text(0.08,0.5, "N", ha='center', va='center',rotation='vertical') 
#
#for index, molec in enumerate(molecList) : 
#    file = 'B_State/%s/fit_hbond/nw_theta2.his'%(molec) 
#    file2 = 'B_State/%s/fit_hbond/nw_theta2.poly'%(molec) 
#    file3 = 'B_State/%s/fit_hbond/nw_geometry.xvg'%(molec) 
#    
#    headlines = 0 
#    with open(file3) as f :
#        lines = f.readlines() 
#        for line in lines : 
#            if line.startswith('#') or line.startswith('@') : 
#                 headlines += 1 
#            else : 
#                 break 
#
#    data = np.genfromtxt(file) 
#    data2 = np.genfromtxt(file2) 
#    data3 = np.genfromtxt(file3,skip_header=headlines) 
#
#    data[:,1] *= len(data3[:,0]) 
#    data2[:,1] *= len(data3[:,0]) 
#
#    ax = axarr[index/figCols,index%figCols]
#
#    ax.plot(data[:,0],data[:,1])
#    ax.plot(data2[:,0],data2[:,1])
#    ax.set_title(file.split('/')[1]) 
#    ax.set_ylim(0,450  ) 
#    ax.set_xlim(100,180) 
#
#fig.savefig('figures/Geometries_theta2_nw_angles.png',format='png') 
#plt.close() 
