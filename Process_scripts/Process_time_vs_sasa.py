import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 
from matplotlib import rc_file

figCols=3
figRows=3

rcFile = 'rc_files/paper.rc'
rc_file(rcFile) 

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

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.1) 
fig.text(0.5,0.04, r"Time (ns)", ha='center', va='center') 
fig.text(0.03,0.5, r"SASA (nm$^2$)", ha='center', va='center',rotation='vertical') 

index=0
for mol in molecList : 
    print index, index/figCols, index%figRows
    ax = axarr[index/figCols,index%figCols]

    datafile = 'B_State/GFP_%s/sasa/cnf_area.xvg'%mol
    print datafile
    try : 
        data1 = np.genfromtxt(datafile,skip_header=27) 
        data1[:,0] = data1[:,0] / 1000 
        equilTime = len(data1) / 5 
        data1 = data1[equilTime:] ##Discard first 10ns as equilibration time
        ax.scatter(data1[:,0],data1[:,2],color=nameToColorKeys[mol],s=0.1) 
    except IOError : 
        print "No file found for %s"%(datafile)  
    except : 
        print "Error importing data from file %s"%(datafile)

    ax.set_title(datafile.split('/')[1].split('_')[1]) 
    ax.set_xlim(10,50) 
    ax.set_ylim(0,2.5)     

    index +=1

fig.savefig('figures/sasa_v_time.png',format='png') 
plt.close() 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15) 
fig.text(0.5,0.04, r"Time (ns)", ha='center', va='center') 
fig.text(0.03,0.5, r"SASA (nm$^2$)", ha='center', va='center',rotation='vertical') 

index=0
for mol in molecList : 
    print index, index/figCols, index%figRows
    ax = axarr[index/figCols,index%figCols]

    datafile = 'B_State/GFP_%s/sasa/aromatic.xvg'%mol
    print datafile
    try : 
        data1 = np.genfromtxt(datafile,skip_header=27) 
        data1[:,0] = data1[:,0] / 1000 
        equilTime = len(data1) / 5 
        data1 = data1[equilTime:] ##Discard first 10ns as equilibration time
        ax.scatter(data1[:,0],data1[:,2],color=nameToColorKeys[mol],s=0.1) 
    except IOError : 
        print "No file found for %s"%(datafile)  
    except : 
        print "Error importing data from file %s"%(datafile)

    ax.set_title(datafile.split('/')[1].split('_')[1]) 
    ax.set_xlim(10,50) 
    ax.set_ylim(0,2.5)     

    index +=1

fig.savefig('figures/sasa_v_time_aromatic.png',format='png') 
plt.close() 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15) 
fig.text(0.5,0.04, r"Time (ns)", ha='center', va='center') 
fig.text(0.03,0.5, r"SASA (nm$^2$)", ha='center', va='center',rotation='vertical') 

index=0
for mol in molecList : 
    print index, index/figCols, index%figRows
    ax = axarr[index/figCols,index%figCols]

    datafile = 'B_State/GFP_%s/sasa/nh_ct_area.xvg'%mol
    print datafile
    try : 
        data1 = np.genfromtxt(datafile,skip_header=27) 
        data1[:,0] = data1[:,0] / 1000 
        equilTime = len(data1) / 5 
        data1 = data1[equilTime:] ##Discard first 10ns as equilibration time
        ax.scatter(data1[:,0],data1[:,2],color=nameToColorKeys[mol],s=0.1) 
    except IOError : 
        print "No file found for %s"%(datafile)  
    except : 
        print "Error importing data from file %s"%(datafile)

    ax.set_title(datafile.split('/')[1].split('_')[1]) 
    ax.set_xlim(10,50) 
    ax.set_ylim(0,0.75)    

    index +=1

fig.savefig('figures/sasa_v_time_nitrile.png',format='png') 
plt.close() 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.1) 
fig.text(0.5,0.04, r"Time (ns)", ha='center', va='center') 
fig.text(0.02,0.5, r"$\left < \rm{SASA} \right > (\rm{nm}^2$)", ha='center', va='center',rotation='vertical') 

figD, axarrD = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
figD.subplots_adjust(wspace=0.1,hspace=0.35,left=0.1) 
figD.text(0.5,0.04, r"Time (ns)", ha='center', va='center') 
figD.text(0.02,0.5, r"$\frac{d}{dt} \left < \rm{SASA} \right > (\rm{nm}^2 \rm{nm}^{-1}$)", ha='center', va='center',rotation='vertical') 

index=0
for mol in molecList : 
    print index, index/figCols, index%figRows
    ax = axarr[index/figCols,index%figCols]
    axD = axarrD[index/figCols,index%figCols]

    datafile = 'B_State/GFP_%s/sasa/cnf_area.xvg'%mol
    try : 
        data1 = np.genfromtxt(datafile,skip_header=27) 
        data1 = data1[2500:,:] 
        data1[:,0] = data1[:,0] / 1000 
    except IOError : 
        print "No file found for %s"%(datafile)  
    except ValueError :
        print "Trying %s again without bottom line."%datafile
        data1 = np.genfromtxt(datafile,skip_header=27,skip_footer=1) 
        data1 = data1[2500:,:] 
        data1[:,0] = data1[:,0] / 1000 
    
    avg, firstD = [],[]
    for i in range(0,len(data1[:,2]) ) : 
        avg.append(np.average(data1[:i,2])) 
    for i in range(len(avg) - 1) : 
        firstD.append(avg[i+1] - avg[i]) 
    firstD = np.array(firstD) 

    ax.plot(data1[:,0],avg,color='b') 
    firstD= firstD
    axD.plot(data1[1:,0],firstD,color='b',linewidth=0.1) 

#    datafile = mol+'/sidechain.xvg'
#    try : 
#        data1 = np.genfromtxt(datafile,skip_header=27) 
#        data1 = data1[2500:,:] 
#        data1[:,0] = data1[:,0] / 1000 
#    except IOError : 
#        print "No file found for %s"%(datafile)  
#    except ValueError :
#        print "Trying %s again without bottom line."%datafile
#        data1 = np.genfromtxt(datafile,skip_header=27,skip_footer=1) 
#        data1 = data1[2500:,:] 
#        data1[:,0] = data1[:,0] / 1000 
#    
#    avg, firstD = [],[]
#    for i in range(0,len(data1[:,2]) ) : 
#        avg.append(np.average(data1[:i,2])) 
#    for i in range(len(avg) - 1) : 
#        firstD.append(avg[i+1] - avg[i]) 
#
#    ax.plot(data1[:,0],avg,color='r') 
#    axD.plot(data1[1:,0],firstD,color='r',linewidth=0.1) 

#    datafile = mol+'/nit_4_atoms_area.xvg'
#    try : 
#        data1 = np.genfromtxt(datafile,skip_header=27) 
#        data1 = data1[2500:,:] 
#        data1[:,0] = data1[:,0] / 1000 
#    except IOError : 
#        print "No file found for %s"%(datafile)  
#    except ValueError :
#        print "Trying %s again without bottom line."%datafile
#        data1 = np.genfromtxt(datafile,skip_header=27,skip_footer=1) 
#        data1 = data1[2500:,:] 
#        data1[:,0] = data1[:,0] / 1000 
#    
#    avg, firstD = [],[]
#    for i in range(0,len(data1[:,2]) ) : 
#        avg.append(np.average(data1[:i,2])) 
#    for i in range(len(avg) - 1) : 
#        firstD.append(avg[i+1] - avg[i]) 
#    firstD = np.array(firstD) 
#
#    ax.plot(data1[:,0],avg,color='g') 
#    firstD= firstD - 0.0005
#    axD.plot(data1[1:,0],firstD,color='g',linewidth=0.1) 

#    datafile = mol+'/nh_ct_area.xvg'
#    try : 
#        data1 = np.genfromtxt(datafile,skip_header=27) 
#        data1 = data1[2500:,:] 
#        data1[:,0] = data1[:,0] / 1000 
#    except IOError : 
#        print "No file found for %s"%(datafile)  
#    except ValueError :
#        print "Trying %s again without bottom line."%datafile
#        data1 = np.genfromtxt(datafile,skip_header=27,skip_footer=1) 
#        data1 = data1[2500:,:] 
#        data1[:,0] = data1[:,0] / 1000 
#    
#    avg, firstD = [],[]
#    for i in range(0,len(data1[:,2]) ) : 
#        avg.append(np.average(data1[:i,2])) 
#    for i in range(len(avg) - 1) : 
#        firstD.append(avg[i+1] - avg[i]) 
#
#    ax.plot(data1[:,0],avg,color='b') 
#    axD.plot(data1[1:,0],firstD,color='b',linewidth=0.1) 



    ax.set_title(datafile.split('/')[1].split('_')[1]) 
    ax.set_xlim(10,50) 
    ax.set_ylim(0,2.5)     

    axD.set_title(datafile.split('/')[1].split('_')[1]) 
    axD.set_xlim(10,50) 
    axD.set_ylim(-.001,.001)    

    index +=1

figD.savefig('figures/sasa_v_time_firstD.png',format='png') 
fig.savefig('figures/sasa_v_time_avg.pdf',format='pdf') 
plt.close() 
