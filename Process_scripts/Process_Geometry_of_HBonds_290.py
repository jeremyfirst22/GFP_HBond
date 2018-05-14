import glob 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm
import os
import matplotlib as mpl 
from scipy.stats import linregress
from matplotlib import rc_file

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

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.10,hspace=0.35,left=0.15,right=0.98,top=0.93,bottom=0.1) 
fig.text(0.5,0.02, r"$d_{\rm{NH}}$ ($\rm{\AA}$)", ha='center', va='center') 
fig.text(0.05,0.5, r"$\theta_1$ (deg)", ha='center', va='center',rotation='vertical') 

xmin,xmax = 1.45, 2.45
ymin,ymax = 99,180 

for index,molec in enumerate(molecList) : 
    file="B_State/GFP_%s/hbond_with_ca_290/geometry.xvg"%molec
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

    ax.set_title(molec) 
    ax.set_ylim([ymin,ymax])
    ax.set_xlim([xmin,xmax]) 

fig.savefig('figures/Geometries_heat_290.png',format='png') 
plt.close() 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.10,hspace=0.35,left=0.15,right=0.98,top=0.93,bottom=0.1) 
fig.text(0.5,0.02, r"$d_{\rm{NH}}$ ($\rm{\AA}$)", ha='center', va='center') 
fig.text(0.05,0.5, r"$\theta_1$ (deg)", ha='center', va='center',rotation='vertical') 

xmin,xmax = 1.45, 2.45
ymin,ymax = 99,180 

for index,molec in enumerate(molecList) : 
    file="B_State/GFP_%s/hbond_with_ca_290/nw_geometry.xvg"%molec
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

    ax.set_title(molec) 
    ax.set_ylim([ymin,ymax])
    ax.set_xlim([xmin,xmax]) 

fig.savefig('figures/Geometries_protein_heat_290.png',format='png') 
plt.close() 

###Plot histogrammed data
#
#
###Theta 1 analyis
fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.10,hspace=0.35,left=0.15,right=0.98,top=0.93,bottom=0.1) 
fig.text(0.5,0.02, r"$\theta_1$", ha='center', va='center') 
fig.text(0.04,0.5, r"Occurances", ha='center', va='center',rotation='vertical') 

figDens, axarrDens = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
figDens.subplots_adjust(wspace=0.10,hspace=0.35,left=0.15,right=0.98,top=0.93,bottom=0.1) 
figDens.text(0.5,0.04, r"$\theta_1$", ha='center', va='center') 
figDens.text(0.04,0.5, r"$\rho_{\rm{dens}}$ Occurances/$\rm{\AA}^3$", ha='center', va='center',rotation='vertical') 

f2, ax2 = plt.subplots(1,1) 
f2.subplots_adjust(left=0.20,right=0.95,top=0.95) 
f2.text(0.5,0.04, r"$\theta_1$", ha='center', va='center') 
f2.text(0.04,0.5, r"$\tilde{\nu}$ (cm$^{-1}$)", ha='center', va='center',rotation='vertical') 

absMax = {} 
with open('Exp_data/sasa_abs_data_290.dat') as f : 
    for line in f : 
        if not line.startswith('#') : 
            key = line.split()[0] 
            value = float(line.split()[2]) 
            absMax[key] = value 


accumAvgTheta = [] 
accumAbsMax = [] 
for index, molec in enumerate(molecList) : 
    fileWhis = 'B_State/GFP_%s/fit_hbond_with_ca_290/theta1.his'%(molec) 
    fileWpoly = 'B_State/GFP_%s/fit_hbond_with_ca_290/theta1.poly'%(molec) 
    fileWgeom = 'B_State/GFP_%s/fit_hbond_with_ca_290/geometry.xvg'%(molec) 
    filePhis = 'B_State/GFP_%s/fit_hbond_with_ca_290/nw_theta1.his'%(molec) 
    filePpoly = 'B_State/GFP_%s/fit_hbond_with_ca_290/nw_theta1.poly'%(molec) 
    filePgeom = 'B_State/GFP_%s/fit_hbond_with_ca_290/nw_geometry.xvg'%(molec) 
    
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
        volumesW[i] = 2*np.pi * r**3 / 3 * (-np.cos((anglesW[i]+binSize/2) * np.pi / 180)  + np.cos((anglesW[i]-binSize/2)* np.pi / 180.) )
    print volumesW
    probsW = probsW / volumesW

    volumesP = np.zeros_like(anglesP) 
    for i in range(len(volumesP)) : 
        if not i == len(volumesP) - 1 : 
            binSize = anglesP[i+1] - anglesP[i]
        else : 
            binSize = anglesP[i] - anglesP[i-1]
        volumesP[i] = 2*np.pi * r**3 / 3 * (-np.cos((anglesP[i]+binSize/2) * np.pi / 180)  + np.cos((anglesP[i]-binSize/2)* np.pi / 180.) )
    print volumesP
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

    if not molec == "Y143X" : 
    #if True : 
        try : 
            ax2.scatter(avgAngleTot, absMax["GFP_"+molec],marker='P',label=molec,color=nameToColorKeys[molec],edgecolor='k',s=100)  
            accumAbsMax.append(absMax["GFP_"+molec]) 
            accumAvgTheta.append(avgAngleTot) 

        except KeyError : 
            print "No key found for %s"%file.split('/')[1]
    else : 
        ax2.scatter(avgAngleTot, absMax["GFP_"+molec],marker='P',label=molec,color=nameToColorKeys[molec],edgecolor='k',s=100,alpha=0.25)  

    ax2.axhline(2227.5,color='k',linestyle='--') 
    ax2.axhline(2235.9,color='b',linestyle='--') 
    
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

    ax.set_title(molec) 
    ax.set_ylim(0,450  ) 
    ax.set_xlim(100,180) 

    axD.set_title(molec) 
    axD.set_ylim(0,3500 ) 
    axD.set_xlim(100,180) 

slope,intercept,r_value,p_value,std_error = linregress(accumAvgTheta, accumAbsMax)
print "r = %f, p = %f"%(r_value,p_value) 

xs = np.linspace(min(accumAvgTheta), max(accumAvgTheta) ) 
ys = slope * xs + intercept
ax2.plot(xs, ys, label="r = %0.3f"%r_value)

ax2.legend(loc=4) 
fig.savefig('figures/Geometries_angles_290.png',format='png') 
figDens.savefig('figures/Geometries_angles_density_290.png',format='png') 
f2.savefig('figures/abs_max_vs_max_theta_290.pdf',format='pdf') 
plt.close() 

