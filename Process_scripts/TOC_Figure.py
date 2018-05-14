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

#rc_file(rcFile) 

if not os.path.isdir('figures') : 
    os.mkdir('figures') 


#

plt.xkcd()

f2, ax2 = plt.subplots(1,1,figsize=(1.5,.875 ) ) 
f2.subplots_adjust(left=0.10,bottom=0.13,right=0.95,top=0.95) 
#f2.text(0.5,0.05, "Angle", ha='center', va='center') 
#f2.text(0.05,0.5, "Frequency",ha='center', va='center',rotation='vertical') 

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

    #ax = axarr[index/figCols,index%figCols]
    #axD = axarrDens[index/figCols,index%figCols]

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

    



#ax2.legend(loc=4,edgecolor='k',framealpha=1.0) 

slope,intercept,r_value,p_value,std_error = linregress(accumAvgTheta, accumAbsMax)

xs = np.linspace(min(accumAvgTheta), max(accumAvgTheta) ) 
ys = slope * xs + intercept
ax2.plot(xs, ys, color='k')

#f2.text(0.23,0.90,r"r = %0.3f"%r_value) 

#fig.savefig('figures/Geometries_angles.png',format='png') 
#figDens.savefig('figures/Geometries_angles_density.png',format='png') 
plt.xticks([]) 
plt.yticks([]) 
f2.savefig('figures/TOC_Graphic_Fig_10.pdf',format='pdf') 
plt.close() 
