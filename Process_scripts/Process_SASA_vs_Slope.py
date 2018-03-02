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

names = np.genfromtxt(inFile,skip_header=1,usecols=0,dtype='str')
slopes= np.genfromtxt(inFile,skip_header=1,usecols=1) 
sasas = np.zeros(len(slopes))
stds = np.zeros(np.size(sasas) ) 

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

for index, name in enumerate(names) : 
    fileName="B_State/%s/sasa/cnf_area.xvg"%(name) 
    if os.path.isfile(fileName) : 
        data = np.genfromtxt(fileName,skip_header=23)
        ##Step size is 4 ps. 2500 frames is 10 ns. Discarded as equiibration time. 
        sasas[index] = np.average(data[equilTime:,2]) 
        stds[index] = np.std(data[equilTime:,2]) 
    else : 
        print "%s file not found! "%fileName 
        continue 


    print name,"\t",sasas[index],"\t",slopes[index],"\t",nameToColorKeys[name]
    plt.errorbar(sasas[index],slopes[index],label=name.split('_', 1)[1],marker='o',xerr=stds[index],color=nameToColorKeys[name],markeredgecolor='k',capsize=3)

index = 0
for name in names : 
    if (name == 'GFP_Y92X' or name == 'GFP_M218X') :  ##Y92X is hbonding to protein M218X is extremely restricted
        sasas = np.delete(sasas,index, None)
        slopes= np.delete(slopes,index, None)
        print "%s deleted"%name 
        continue
    index+=1

slope,intercept,r_value,p_value,std_error = linregress(sasas,slopes) 
print "slope = ",slope,"\tr_value = ",r_value 


x = np.arange(0,3.0,0.01)
plt.ylim([-0.065,-0.01]) 
plt.plot(x,x*slope+intercept,label="r = %0.3f"%r_value)

box = ax.get_position() 
ax.set_position([box.x0, box.y0, box.width *0.85, box.height]) 

plt.legend(bbox_to_anchor=(1.40, 0.95),numpoints=1,edgecolor='k')
plt.xlabel(r"SASA (nm$^2$)") 
plt.ylabel(r"FTLS (cm$^{-1}$ / K)") 


plt.savefig("%s/sasa_vs_slope.pdf"%(saveDir),format='pdf') 
plt.close() 

