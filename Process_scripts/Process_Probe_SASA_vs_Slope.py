import numpy as np 
import matplotlib.pyplot as plt 
import glob as glob 
import os 
from os import sys
from scipy.stats import linregress
from matplotlib import rc_file

inFile='Exp_data/sasa_abs_data.dat'  
rcFile='paper.rc'

saveDir='figures/sasa_vs_slopes' 

if not os.path.isdir(saveDir) : 
    os.mkdir(saveDir)

## Plot for both states
#rc_file("rc_files/%s"%rcFile) 
fig, ax = plt.subplots(1,1) 

names = np.genfromtxt(inFile,skip_header=1,usecols=0,dtype='str')
slopes= np.genfromtxt(inFile,skip_header=1,usecols=1) 
sasas = np.zeros(len(slopes))
stds = np.zeros(np.size(sasas) ) 

for index, name in enumerate(names) : 
    fileName="B_State/%s/sasa/nh_ct_area.xvg"%(name) 
    if os.path.isfile(fileName) : 
        try :
            data = np.genfromtxt(fileName,skip_header=23)
        except ValueError :
            print "Trying %s withou last line"%fileName
            data = np.genfromtxt(fileName,skip_header=23,skip_footer=1) 
        ##Step size is 4 ps. 2500 frames is 10 ns. Discarded as equiibration time. 
        sasas[index] = np.average(data[2500:,2]) 
        stds[index] = np.std(data[2500:,2]) 
    else : 
        print "%s file not found! "%fileName 
        continue 


    print name,"\t",sasas[index],"\t",slopes[index]
    plt.errorbar(sasas[index],slopes[index],label=name.split('_', 1)[1],marker='o',xerr=stds[index])

index = 0
for name in names : 
     if (name == 'GFP_Y92X' or name == 'GFP_M218X' or name == 'GFP_D190X' ) : 
    #if (name == 'GFP_Y92X') : 
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
ax.set_position([box.x0, box.y0, box.width *0.75, box.height]) 

plt.legend(bbox_to_anchor=(1.5, 0.95),numpoints=1,fontsize='x-small')
plt.xlabel(r"SASA (nm$^2$)") 
plt.ylabel(r"$\frac{\partial \tilde{\nu}}{\partial T}$ (cm$^{-1}$ / K)") 


plt.savefig("%s/probe_sasa_vs_slope.pdf"%saveDir,format='pdf') 
plt.close() 

