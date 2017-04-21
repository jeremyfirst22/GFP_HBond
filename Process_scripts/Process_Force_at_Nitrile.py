import glob 
import numpy as np 
import matplotlib.pyplot as plt
from sys import exit
import os 
from scipy.stats import linregress

from matplotlib import rc_file 

force_dir='force_calc'
save_dir='figures/force_calc'
experimental_data='sasa_abs_data.dat'

if not os.path.isdir(save_dir) : 
    os.mkdir(save_dir) 

## Gather experimental data into dictionaries with molecule names
mutToExp = {} 
with open(experimental_data,'r') as f : 
    for line in f : 
        if line[0] == '#' or line[0] == '@' : continue 
        try : 
            key  =line.split()[0]
            value=line.split()[2]
        except : 
            print "WARNING: Key value pair import failed for line: %s"%line
            continue 
        mutToExp[key] = value 

mutExppKa={
    'GFP_Y92X':[7.5,0.5],
    'GFP_F114X':[7.5,0.5],
    'GFP_Y143X':[7.5,0.5],
    'GFP_F145X':[7.5,0.5],
    'GFP_Y149X':[7.5,0.5],
    'GFP_F165X':[7.5,0.5],
    'GFP_N212X':[7.5,0.5],
    'GFP_M218X':[7.5,0.5]
} 

for force_type in ['external_field','solvent_rxn_field','total_field'] : 
    ##Fine all force files 
    force_files = glob.glob('[AB]_State/*/%s/*.%s.projected.xvg'%(force_dir,force_type)) 
    
    molList= []
    for item in force_files : 
        state, molec,throwAway,field = item.split('/') 
        if not os.path.isfile("%s/%s/%s/%s.%s.projected.xvg"%(state,molec,force_dir,molec,force_type)) : 
            print "%s\tNot found, removing!"%molec
            force_files.remove(item) 
        else : 
            molList.append([state, molec]) 
    
    molList = sorted(molList ) 
    numMols = len(molList) /2
    data = [] 
    for molec in range(numMols ) : 
        if not molList[molec][1] == molList[molec+numMols][1] : 
            print "Logical fallicy! Your A & B states do not align. Cowardly exitting." 
            exit() 
        try : 
            fieldA = np.genfromtxt("%s/%s/%s/%s.%s.projected.xvg"%(molList[molec][0],molList[molec][1],force_dir,molList[molec][1],force_type),comments='#' ) 
            fieldB = np.genfromtxt("%s/%s/%s/%s.%s.projected.xvg"%(molList[molec+numMols][0],molList[molec+numMols][1],force_dir,molList[molec+numMols][1],force_type),comments='#' ) 
        except : 
            print "ERROR: %s failed!"%molList[molec]
            continue 
        try : 
            print mutToExp[molList[molec][1]]
        except : 
            print molList[molec][1]
            continue
        
        ##Step size is 4 ps. 2500 frames is 10 ns. Discarded as equilibration time.
        solventFieldA = np.average(fieldA[2500:]) 
        solventFieldB = np.average(fieldB[2500:]) 
        stdA = np.std(fieldA[2500:])
        stdB = np.std(fieldB[2500:])
        
        data.append([molList[molec][1],solventFieldA, solventFieldB, stdA, stdB]) 
    
    for molec in data : 
        try : 
            pKa = float(mutExppKa[molec[0]][0]) 
        except (KeyError, ValueError) : 
            print "%s pKa not found!"%molec[0]
            data.remove(molec) 
            print molec[0], " removed!"
            continue
        #print molec[0],pKa
        ## From Henderson-Hasselbalch eqn. 
        ratio = 10**+(7.4 - pKa) 
        wA = 1/(1+ratio) 
        wB = ratio / (ratio +1) 
        weightedF = molec[2] * wB + molec[1]*wA 
        weightedSTD = np.sqrt(molec[3]**2/wA + molec[4]**2/wB ) 
        molec.append(weightedF) 
        molec.append(weightedSTD)
    
    
    ## Print calculated data to a file
    force_file="%s/%s.dat"%(save_dir,force_type) 
    with open(force_file,'w') as f: 
        f.write("Name      SolventFieldA   SolventFieldB    WeightedField   WeightedSTD\n") 
    with open(force_file,'a') as f : 
        for item in data : 
            f.write("%s\t%f\t%f\t%f\t%f\n"%(item[0],item[1],item[2],item[5],item[6]) ) 
    
    ## Read in data from file and plot 
    names= np.genfromtxt(force_file,skip_header=1,usecols=0,dtype='str') 
    dataA= np.genfromtxt(force_file,skip_header=1,usecols=1) 
    dataB= np.genfromtxt(force_file,skip_header=1,usecols=2) 
    dataW= np.genfromtxt(force_file,skip_header=1,usecols=3) 
    dataS= np.genfromtxt(force_file,skip_header=1,usecols=4) 
    
    fig1, ax1 = plt.subplots(1,1) 
    for index, name in enumerate(names) : 
        print name,"\t",mutToExp[name],"\t",dataW[index],"\t",dataA[index],"\t",dataB[index]
        ax1.errorbar(mutToExp[name],dataW[index],label=name,marker='o')
        ax1.scatter(mutToExp[name],dataA[index],marker='D',color='darkgray')
        ax1.scatter(mutToExp[name],dataB[index],marker='^',color='darkgray')

    #minX, maxX = np.min(dataW), np.max(dataW) 
    #x = np.arange(minX,maxX,0.01) 
    
    fig1.savefig("%s/%s.pdf"%(save_dir,force_type),format='pdf') 
    plt.close() 


