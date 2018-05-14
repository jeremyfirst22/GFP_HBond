import glob 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm
from scipy.stats import linregress
import os
from sys import exit
from matplotlib import rc_file

rcFile='rc_files/paper.rc'
rc_file(rcFile) 

figCols=3
figRows=3

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

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

FTLS = {}  
with open('Exp_data/sasa_abs_data_290.dat') as f : 
    for line in f : 
        if not line.startswith('#') : 
            key = line.split()[0] 
            value = float(line.split()[1]) 
            FTLS[key] = value 

molecKList = {} 
for temp in [290,300,310] : 
    for index,molec in enumerate(molecList) :
        if not temp == 300 : file = 'B_State/GFP_%s/hbond_with_ca_%i/hb_count.xvg'%(molec,temp) 
        else : file = 'B_State/GFP_%s/hbond_with_ca/hb_count.xvg'%molec  
        data = np.genfromtxt(file,skip_header=25) 
        hbonds = data[:,6] 
         
        equilK = (np.sum(data[:,6]) - data[0,6]) / np.sum(data[:,6])*100
    
        molec = file.split('/')[1].split('_')[1]
        print molec, temp, equilK

#        ax = axarr[index/figCols,index%figCols]

        if not molec in molecKList : 
            molecKList[molec] = [[temp, equilK]] 
        else : 
            molecKList[molec].append([temp, equilK]) 
        
    
print molecKList

f1, ax1 = plt.subplots(1,1)
f1.subplots_adjust(left=0.20,right=0.95,top=0.95)
f1.text(0.5,0.04, r"1/T",ha='center', va='center')
f1.text(0.04,0.5, r"ln(K)",ha='center', va='center',rotation='vertical')

FTLSSlopes = {} 
for molec in molecKList : 
    data = np.array(molecKList[molec]) 
    x = data[:,0]
    y = data[:,1]
    ax1.scatter(1/x, np.log(y),color=nameToColorKeys[molec],label=molec) 

    slope, intercept, r_value, p_value, std_error = linregress(1/x, np.log(y)) 
    xs = np.linspace(min(1/x), max(1/x) ) 
    ys = slope * xs + intercept 
    ax1.plot(xs,ys,color=nameToColorKeys[molec]) 
    
    FTLSSlopes[molec] = [slope, float(FTLS["GFP_"+molec]), std_error]
    

ax1.legend(loc=4) 
ax1.set_xlim(1/320., 1/275.) 
f1.savefig('figures/Vant_hoff.pdf',format='pdf') 

f2, ax2 = plt.subplots(1,1)
f2.subplots_adjust(left=0.20,right=0.95,top=0.95)
f2.text(0.5,0.04, r"$-\frac{\Delta H}{R}$",ha='center', va='center')
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
ax2.legend(loc=4) 

f2.savefig('figures/Vant_Hoff_Slopes_vs_FTLS.png',format='png')

