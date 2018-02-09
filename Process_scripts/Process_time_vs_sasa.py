import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 

figCols=3
figRows=3

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0) 
fig.text(0.5,0.04, r"Time (ns)", ha='center', va='center') 
fig.text(0.04,0.5, r"SASA (nm$^2$)", ha='center', va='center',rotation='vertical') 

molList = glob.glob('B_State/*/sasa') 
index=0
for mol in molList : 
    print index, index/figCols, index%figRows
    ax = axarr[index/figCols,index%figCols]

    datafile = mol+'/cnf_area.xvg'
    print datafile
    try : 
        data1 = np.genfromtxt(datafile,skip_header=27) 
        data1[:,0] = data1[:,0] / 1000 
        ax.scatter(data1[:,0],data1[:,2],color='k',s=0.1) 
    except IOError : 
        print "No file found for %s"%(datafile)  
    except ValueError :
        print "Trying %s again without bottom line."%datafile
        data1 = np.genfromtxt(datafile,skip_header=27,skip_footer=1) 
        data1[:,0] = data1[:,0] / 1000 
        ax.scatter(data1[:,0],data1[:,2],color='k',s=0.1,label="Whole residue") 

    #datafile = mol+'/sidechain.xvg'
    #print datafile
    #try : 
    #    data1 = np.genfromtxt(datafile,skip_header=27) 
    #    data1[:,0] = data1[:,0] / 1000 
    #    ax.scatter(data1[:,0],data1[:,2],color='r',s=0.1,label="Sidechain only") 
    #except IOError : 
    #    print "No file found for %s"%(datafile)  
    #except ValueError :
    #    print "Trying %s again without bottom line."%datafile
    #    data1 = np.genfromtxt(datafile,skip_header=27,skip_footer=1) 
    #    data1[:,0] = data1[:,0] / 1000 
    #    ax.scatter(data1[:,0],data1[:,2],color='r',s=0.1) 

    datafile = mol+'/nit_4_atoms_area.xvg'
    try : 
        data2 = np.genfromtxt(datafile,skip_header=27) 
        data2[:,0] = data2[:,0] / 1000 
        ax.scatter(data2[:,0],data2[:,2],color='g',s=0.1,label="5 atom area") 
    except IOError : 
        print "No file found for %s"%(datafile)  
    except ValueError :
        print "Trying %s again without bottom line."%datafile
        data2 = np.genfromtxt(datafile,skip_header=27,skip_footer=1) 
        data2[:,0] = data2[:,0] / 1000 
        ax.scatter(data2[:,0],data2[:,2],color='g',s=0.1) 

    #datafile = mol+'/nh_ct_area.xvg'
    #try : 
    #    data3 = np.genfromtxt(datafile,skip_header=27) 
    #    data3[:,0] = data3[:,0] / 1000 
    #    ax.scatter(data3[:,0],data3[:,2],color='b',s=0.1,label="C N atoms only" ) 
    #except IOError : 
    #    print "No file found for %s"%(datafile)  
    #except ValueError :
    #    print "Trying %s again without bottom line."%datafile
    #    data3 = np.genfromtxt(datafile,skip_header=27,skip_footer=1) 
    #    data3[:,0] = data3[:,0] / 1000 
    #    ax.scatter(data3[:,0],data3[:,2],color='b',s=0.1,label="C N atoms only") 


    ax.set_title(datafile.split('/')[1].split('_')[1]) 
    ax.set_xlim(0,50) 
    ax.set_ylim(0,3)     

    index +=1

fig.savefig('figures/sasa_v_time.png',format='png') 
plt.close() 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0) 
fig.text(0.5,0.04, r"Time (ns)", ha='center', va='center') 
fig.text(0.02,0.5, r"$\left < \rm{SASA} \right > (\rm{nm}^2$)", ha='center', va='center',rotation='vertical') 

figD, axarrD = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
figD.subplots_adjust(wspace=0) 
figD.text(0.5,0.04, r"Time (ns)", ha='center', va='center') 
figD.text(0.02,0.5, r"$\frac{d}{dt} \left < \rm{SASA} \right > (\rm{nm}^2 \rm{nm}^{-1}$)", ha='center', va='center',rotation='vertical') 

molList = glob.glob('B_State/*/sasa') 
index=0
for mol in molList : 
    print index, index/figCols, index%figRows
    ax = axarr[index/figCols,index%figCols]
    axD = axarrD[index/figCols,index%figCols]

    datafile = mol+'/cnf_area.xvg'
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

    ax.plot(data1[:,0],avg,color='k') 
    firstD= firstD + 0.0005
    axD.plot(data1[1:,0],firstD,color='k',linewidth=0.1) 

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

    datafile = mol+'/nit_4_atoms_area.xvg'
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

    ax.plot(data1[:,0],avg,color='g') 
    firstD= firstD - 0.0005
    axD.plot(data1[1:,0],firstD,color='g',linewidth=0.1) 

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
