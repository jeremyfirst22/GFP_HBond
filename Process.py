import glob as glob 
import numpy as np 
import matplotlib.pyplot as plt 

datafiles = glob.glob('GFP_*X/HBond/cnf_num.xvg') 
print datafiles 

for file in datafiles : 
    data = np.genfromtxt(file,skip_header=23) 
    print "%10s\t%5f\t%5f"%(file.split('/')[0], np.average(data[:,1]), np.std(data[:,1]) ) 
