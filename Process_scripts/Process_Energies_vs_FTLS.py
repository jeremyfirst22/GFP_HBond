import numpy as np 
import matplotlib.pyplot as plt 
import glob as glob 
import os 
from os import sys
from scipy.stats import linregress
from matplotlib import rc_file

equilTime=10   #ns 
equilTime=equilTime*1000 / 4 ##frames

inFile='Exp_data/sasa_abs_data.dat'  
#rcFile='paper.rc'

saveDir='figures'

if not os.path.isdir(saveDir) : 
    os.mkdir(saveDir)

#rc_file("rc_files/%s"%rcFile) 
fig, ax = plt.subplots(1,1) 

names = np.genfromtxt(inFile,skip_header=1,usecols=0,dtype='str')
slopes= np.genfromtxt(inFile,skip_header=1,usecols=1) 

#kvalues = np.genfromtxt('water_chosen_decay.xvg') 
kvalues = {} 
cvalues = {} 
with open('fits/protein_chosen_decay.xvg') as f : 
    for line in f : 
        key = line.split()[0]
        k = line.split()[2]
        c = line.split()[1]
        kvalues[key] = k 
        cvalues[key] = c 

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


energies = [] 
for index, name in enumerate(names) : 
    fileName="B_State/%s/sasa/nh_ct_area.xvg"%(name) 
    if os.path.isfile(fileName) : 
        data = np.genfromtxt(fileName,skip_header=23)
        ##Step size is 4 ps. 2500 frames is 10 ns. Discarded as equiibration time. 
        sasas[index] = np.average(data[equilTime:,2]) 
        stds[index] = np.std(data[equilTime:,2]) 
    else : 
        print "%s file not found! "%fileName 
        continue 


    #kA = float(kvalues[name]) / float(sasas[index] ) / 100 
    kA = float(kvalues[name]) / float(cvalues[name] ) 
    energy = -1 *  np.log(kA ) * 8.314 * 298 
    energies.append(energy) 

    plt.errorbar(energy ,slopes[index],label=name.split('_', 1)[1],marker='o',xerr=stds[index],color=nameToColorKeys[name],markeredgecolor='k',capsize=3)

energies = np.array(energies) 
#index = 0
#for name in names : 
#    if (name == 'GFP_Y92X' or name == 'GFP_Y143X') :  ##Y92X is hbonding to protein M218X is extremely restricted
#        sasas = np.delete(sasas,index, None)
#        energies = np.delete(energies,index, None)
#        slopes= np.delete(slopes,index, None)
#        print "%s deleted"%name 
#        continue
#    index+=1

slope,intercept,r_value,p_value,std_error = linregress(energies,slopes) 
print "slope = ",slope,"\tr_value = ",r_value 


x = np.linspace(energies.min() ,energies.max() ,100 )
#plt.ylim([-0.065,-0.01]) 
#plt.xlim([18000,24000]) 
plt.plot(x,x*slope+intercept,label="r = %0.3f"%r_value)

box = ax.get_position() 
ax.set_position([box.x0, box.y0, box.width *0.85, box.height]) 

plt.legend(bbox_to_anchor=(1.25, 0.95),numpoints=1,fontsize='x-small',edgecolor='k')
plt.xlabel(r"$E_a$ (units?)") 
plt.ylabel(r"FTLS (cm$^{-1}$ / K)") 


plt.savefig("%s/Energies.pdf"%(saveDir),format='pdf') 
plt.close() 

