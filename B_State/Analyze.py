import numpy as np 
import glob as glob 
import matplotlib.pyplot as plt 

datafiles = glob.glob('*/rmsd/*_rms.xvg') 

for file in datafiles : 
    data = np.genfromtxt(file,skip_header=16) 
    plt.plot(data[:,0],data[:,1],label=file.split('/')[0]) 

plt.legend(loc=4) 
plt.show() 
plt.close() 

datafiles = glob.glob('*/gyrate/*_gyrate.xvg') 

for file in datafiles : 
    data = np.genfromtxt(file,skip_header=25) 
    plt.plot(data[:,0],data[:,1],label=file.split('/')[0]) 

plt.legend() 
plt.show() 
plt.close() 
