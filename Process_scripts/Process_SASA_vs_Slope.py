import numpy as np 
import matplotlib.pyplot as plt 
import glob as glob 
import os 
from os import sys
from scipy.stats import linregress

inFile='Exp_data/sasa_abs_data.dat'  

saveDir='figures/sasa_vs_slopes' 

if not os.path.isdir(saveDir) : 
    os.mkdir(saveDir)

## Plot for both states
for state in ['A','B'] : 

    names = np.genfromtxt(inFile,skip_header=1,usecols=0,dtype='str')
    slopes= np.genfromtxt(inFile,skip_header=1,usecols=1) 
    sasas = np.zeros(len(slopes))
    stds = np.zeros(np.size(sasas) ) 

    for index, name in enumerate(names) : 
        fileName="%s_State/%s/sasa/%s.cnf.xvg"%(state,name,name)
        if os.path.isfile(fileName) : 
            data = np.genfromtxt(fileName,skip_header=23)
            ##Step size is 4 ps. 2500 frames is 10 ns. Discarded as equiibration time. 
            sasas[index] = np.average(data[2500:,2]) 
            stds[index] = np.std(data[2500:,2]) 
        else : 
            print "%s file not found! "%fileName 
            sys.exit() 
    
    
        print name,"\t",sasas[index],"\t",slopes[index]
        plt.errorbar(sasas[index],slopes[index],label=name,marker='o',xerr=stds[index])
    
    index = 0
    for name in names : 
        if (name == 'GFP_Y92X' or name == 'GFP_M218X') : 
            sasas = np.delete(sasas,index, None)
            slopes= np.delete(slopes,index, None)
            print "%s deleted"%name 
            continue
        index+=1
    
    slope,intercept,r_value,p_value,std_error = linregress(sasas,slopes) 
    print "slope = ",slope,"\tr_value = ",r_value 
    
    x = np.arange(0,3.0,0.01)
    plt.ylim([-0.065,-0.01]) 
    plt.plot(x,x*slope+intercept,label="r = %f"%r_value)

    plt.legend(loc=1,numpoints=1,fontsize='x-small')
    plt.xlabel(r"SASA (nm$^2$)") 
    plt.ylabel(r"$\frac{\partial \tilde{\nu}}{\partial T}$ (cm$^{-1}$ / K)") 


    plt.savefig("%s/%s_sasa_vs_slope.pdf"%(saveDir,state),format='pdf') 
    plt.close() 

