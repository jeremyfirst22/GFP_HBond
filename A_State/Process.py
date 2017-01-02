import glob as glob 
import numpy as np 
import matplotlib.pyplot as plt 
from scipy.stats import linregress

datafiles = glob.glob('GFP_*X/HBond_nit/*.frame_hb.xvg') 
#print datafiles 


exp_slopes={
'GFP_Y92X' :   -0.060,
'GFP_F114X':   -0.028,
'GFP_Y143X':   -0.033,
'GFP_F145X':   -0.020,
'GFP_Y149X':   -0.023,
'GFP_F165X':   -0.016,
'GFP_N212X':   -0.053,
'GFP_M218X':   -0.035,
}

exper, avgs = [], [] 
for file in datafiles : 
    molec=file.split('/')[0]
    print molec
    data = np.genfromtxt(file,skip_header=23) 
    avg = np.average(data[:,2]) 
    std = np.std(data[:,2]) 
    print "%10s\t%5f\t%5f"%(molec,avg,std) 
    plt.errorbar(exp_slopes[molec],avg,yerr=std,marker='o',label=molec) 
    exper.append(exp_slopes[molec]) 
    avgs.append(avg) 

x = np.arange(-0.060,-0.016,0.01) 
slope, intercept, r_value, p_value, std_err = linregress(exper,avgs) 
plt.plot(x,x*slope+intercept,label="r = %.2f"%r_value) 
plt.legend(loc=4) 


plt.show() 
