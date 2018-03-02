import numpy as np 
import matplotlib.pyplot as plt 
import glob as glob 
import os 
from os import sys
from scipy.stats import linregress
from matplotlib import rc_file
import matplotlib as mpl 

equilTime=10   #ns 
equilTime=equilTime*1000 / 4 ##frames

inFile='Exp_data/sasa_abs_data.dat'  
rcFile='rc_files/paper.rc'

saveDir='figures/sasa_vs_slopes' 

if not os.path.isdir(saveDir) : 
    os.mkdir(saveDir)

rc_file(rcFile) 
mpl.rcParams['figure.figsize'] = 5,4

fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(left=0.15,right=0.85) 

nameToColorKeys = {
        "GFP_F145X":'#A93226',  
        "GFP_M218X":'r',  
        "GFP_F165X":'#E67E22',  
        "GFP_F114X":'y',  
        "GFP_N198X":'#2ECC71',  
        "GFP_Y143X":'g',  
        "GFP_N212X":'c',  
        "GFP_D117X":'b',  
        "GFP_Y92X":'m'  
        }

nameToExp = {} 
with open('Exp_data/FWHM.dat') as f : 
    for line in f : 
        if not line.startswith('#') : 
            key, value = line.split() 
            nameToExp[key] = value 

datafiles = glob.glob('B_State/*/fit_hbond_with_ca/theta1.his') 

fwhmAccum, stdevAccum = [],[]
for dFile in datafiles : 
    molec = dFile.split('/')[1]
    fwhm = float(nameToExp[molec]) 

    pFile = 'B_State/%s/fit_hbond_with_ca/nw_theta1.his'%(molec)
    if molec == "GFP_Y92X" or molec == "GFP_Y143X" or molec == "GFP_F165X" : 
        continue 
        data = np.genfromtxt(pFile) 
        alpha = 0.25
    elif molec == "GFP_F165X" : 
        continue 
    else : 
        data = np.genfromtxt(dFile) 
        alpha = 1 
    angles = data[:,0]
    probs = data[:,1]

    avg = np.average(angles,weights=probs) 
    variance = np.average((angles - avg)**2,weights=probs) 
    stdev = np.sqrt(variance)

    fwhmAccum.append(fwhm) 
    stdevAccum.append(stdev)

    print "%6s %5.1f %3.1f"%(molec, fwhm, stdev) 

    plt.scatter(fwhm,stdev,color=nameToColorKeys[molec],label=molec.split('_')[1],alpha=alpha,edgecolor='k')

slope,intercept,r_value,pvalue,stderr = linregress(fwhmAccum, stdevAccum) 
xs = np.linspace(min(fwhmAccum),max(fwhmAccum),100) 
ys = slope * xs + intercept 

plt.plot(xs,ys,label="r = %.2f"%r_value) 

plt.legend(loc=4,fontsize='xx-small')
plt.savefig('figures/fwhm_vs_stdev.pdf',format='pdf') 


