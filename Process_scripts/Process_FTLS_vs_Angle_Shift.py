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

rc_file(rcFile) 

if not os.path.isdir('figures') : 
    os.mkdir('figures') 


###Plot histogrammed data
#
#
###Theta 1 analyis


FTLS = {} 
with open('Exp_data/sasa_abs_data_290.dat') as f : 
    for line in f : 
        if not line.startswith('#') : 
            key = line.split()[0] 
            value = float(line.split()[1]) 
            FTLS[key] = value 


accumAvgTheta = [] 
accumAbsMax = [] 
avgAngles = {} 
for index, molec in enumerate(molecList) : 
    for temp in [290, 300, 310] : 
        if temp == 300 : 
            fileWhis = 'B_State/GFP_%s/fit_hbond_with_ca/theta1.his'%(molec) 
            fileWpoly = 'B_State/GFP_%s/fit_hbond_with_ca/theta1.poly'%(molec) 
            fileWgeom = 'B_State/GFP_%s/fit_hbond_with_ca/geometry.xvg'%(molec) 
            filePhis = 'B_State/GFP_%s/fit_hbond_with_ca/nw_theta1.his'%(molec) 
            filePpoly = 'B_State/GFP_%s/fit_hbond_with_ca/nw_theta1.poly'%(molec) 
            filePgeom = 'B_State/GFP_%s/fit_hbond_with_ca/nw_geometry.xvg'%(molec) 
        else :  
            fileWhis = 'B_State/GFP_%s/fit_hbond_with_ca_%s/theta1.his'%(molec,str(temp)) 
            fileWpoly = 'B_State/GFP_%s/fit_hbond_with_ca_%s/theta1.poly'%(molec,str(temp)) 
            fileWgeom = 'B_State/GFP_%s/fit_hbond_with_ca_%s/geometry.xvg'%(molec,str(temp)) 
            filePhis = 'B_State/GFP_%s/fit_hbond_with_ca_%s/nw_theta1.his'%(molec,str(temp)) 
            filePpoly = 'B_State/GFP_%s/fit_hbond_with_ca_%s/nw_theta1.poly'%(molec,str(temp)) 
            filePgeom = 'B_State/GFP_%s/fit_hbond_with_ca_%s/nw_geometry.xvg'%(molec,str(temp)) 

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

        dataWhis[:,1] *= len(dataWgeom[:,0]) 
        dataWpoly[:,1] *= len(dataWgeom[:,0]) 

        anglesW = dataWpoly[:,0] 
        probsW = dataWpoly[:,1]

        if not molec == "N212X" : 
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

            dataPhis[:,1] *= len(dataPgeom[:,0]) 
            dataPpoly[:,1] *= len(dataPgeom[:,0]) 

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
            volumesW[i] = np.pi * r**3 / 3 * (-np.cos((anglesW[i]+binSize/2) * np.pi / 180)  + np.cos((anglesW[i]-binSize/2)* np.pi / 180.) )
        #print volumesW
        probsW = probsW / volumesW

        if not molec == "N212X" : 
            volumesP = np.zeros_like(anglesP) 
            for i in range(len(volumesP)) : 
                if not i == len(volumesP) - 1 : 
                    binSize = anglesP[i+1] - anglesP[i]
                else : 
                    binSize = anglesP[i] - anglesP[i-1]
                volumesP[i] = np.pi * r**3 / 3 * (-np.cos((anglesP[i]+binSize/2) * np.pi / 180)  + np.cos((anglesP[i]-binSize/2)* np.pi / 180.) )
            #print volumesP
            probsP = probsP / volumesP

        avgAngleW = 0
        avgAngleP = 0
        for i in range(len(anglesW)) : 
            avgAngleW = avgAngleW + probsW[i] * anglesW[i]

        if not molec == "N212X" : 
            for i in range(len(anglesP)) : 
                avgAngleP = avgAngleP + probsP[i] * anglesP[i]
            avgAngleTot = (avgAngleW + avgAngleP) / (sum(probsP)+sum(probsW))
        else : avgAngleTot = avgAngleW / sum(probsW) 

        avgAngleP = avgAngleP / sum(probsP) 
        avgAngleW = avgAngleW / sum(probsW) 

        if not molec in avgAngles : 
            avgAngles[molec] = [[temp, avgAngleTot]] 
        else : 
            avgAngles[molec].append([temp,avgAngleTot]) 


#print avgAngles 

f1, ax1 = plt.subplots(1,1) 
f1.subplots_adjust(left=0.20,right=0.95,top=0.95) 
f1.text(0.5,0.04, r"$\theta_1$", ha='center', va='center') 
f1.text(0.04,0.5, r"$\tilde{\nu}$ (cm$^{-1}$)", ha='center', va='center',rotation='vertical') 

FTLSSlopes = {} 
for molec in avgAngles : 
    data = avgAngles[molec] 
    data = np.array(data) 

    slope,intercept,r_value,p_value,std_error = linregress(data[:,0],data[:,1]) 
    xs = np.linspace(290,310)
    ys = slope* xs+ intercept

    ax1.plot(xs,ys,color=nameToColorKeys[molec]) 
    ax1.scatter(data[:,0],data[:,1],color=nameToColorKeys[molec])  

    FTLSSlopes[molec] = [slope, float(FTLS["GFP_"+molec]), std_error ] 

f1.savefig('figures/HBond_Angle_vs_Temp.png',format='png')  

f2, ax2 = plt.subplots(1,1) 
f2.subplots_adjust(left=0.20,right=0.95,top=0.95) 
f2.text(0.5,0.04, r"$\frac{d}{dT} <\theta_1>$",ha='center', va='center') 
f2.text(0.04,0.5, r"FTLS",ha='center', va='center',rotation='vertical') 

print FTLSSlopes
xAccum, yAccum = [], [] 
for molec in FTLSSlopes : 
    x,y,x_err = FTLSSlopes[molec]
    ax2.errorbar(x,y,xerr=x_err,color=nameToColorKeys[molec]) 
    ax2.scatter(x,y,color=nameToColorKeys[molec],label=molec) 
    xAccum.append(x)
    yAccum.append(y) 

slope, intercept, r_value,p_value, std_error = linregress(xAccum,yAccum)
x = np.linspace(min(xAccum), max(xAccum)) 
y = slope * x + intercept 
ax2.plot(x,y, color='k', label='r = %f'%r_value) 

ax2.legend(loc=4,fontsize='xx-small') 

f2.savefig('figures/FTLS_vs_HBond_angle_shift.png',format='png') 


