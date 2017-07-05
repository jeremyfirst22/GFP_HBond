import numpy as np 
import matplotlib.pyplot as plt 
import glob as glob 
import os 
from os import sys
from scipy.stats import linregress
from matplotlib import rc_file

inFile='Exp_data/sasa_abs_data.dat'  
rcFile='paper.rc'

saveDir='figures'

if not os.path.isdir(saveDir) : 
    os.mkdir(saveDir)

## Plot for both states
for state in ['A','B'] : 
    rc_file("rc_files/%s"%rcFile) 
    fig, ax = plt.subplots(1,1) 

    names = np.genfromtxt(inFile,skip_header=1,usecols=0,dtype='str')
    slopes= np.genfromtxt(inFile,skip_header=1,usecols=1) 
    sasas = np.zeros(len(slopes))
    stds = np.zeros(np.size(sasas) ) 

    for index, name in enumerate(names) : 
        fileName="%s_State/%s/hbond/rates.txt"%(state,name) 

        with open(fileName) as f : 
            for line in f.readlines() :
                if 'Forward' in line : 
                    rate = line.split()[2] 
        print rate
        rate = float(rate)
    
        plt.errorbar(rate,slopes[index],label=name.split('_', 1)[1],marker='o',xerr=stds[index])
    
    index = 0

    #for name in names : 
#   #     if (name == 'GFP_Y92X' or name == 'GFP_M218X' or name == 'GFP_D190X' ) : 
    #    if (name == 'GFP_Y92X' ) : 
    ##        sasas = np.delete(sasas,index, None)
    #        slopes= np.delete(slopes,index, None)
    #        print "%s deleted"%name 
    #        continue
    #    index+=1
    #
    #slope,intercept,r_value,p_value,std_error = linregress(sasas,slopes) 
    #print "slope = ",slope,"\tr_value = ",r_value 
    
    
    #x = np.arange(0,3.0,0.01)
    #plt.ylim([-0.065,-0.01]) 
    #plt.plot(x,x*slope+intercept,label="r = %0.3f"%r_value)

    box = ax.get_position() 
    ax.set_position([box.x0, box.y0, box.width *0.75, box.height]) 

    plt.legend(bbox_to_anchor=(1.5, 0.95),numpoints=1,fontsize='x-small')
    plt.xlabel(r"Lifetime (ps)") 
    plt.semilogx() 
    plt.ylabel(r"$\frac{\partial \tilde{\nu}}{\partial T}$ (cm$^{-1}$ / K)") 


    plt.savefig("%s/Lifetime_vs_slope_%s.pdf"%(saveDir,state),format='pdf') 
    plt.close() 

