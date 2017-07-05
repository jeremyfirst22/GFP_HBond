import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from os.path import basename
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 

figCols=4
figRows=3

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

molecList=['GFP_Y92X', 'GFP_F114X', 'GFP_Y143X', 'GFP_F145X', 'GFP_Y149X', 'GFP_F165X', 'GFP_N212X', 'GFP_M218X', 'GFP_D117X', 'GFP_D190X', 'GFP_N198X'] 

for state in 'A_State','B_State' : 
    fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
    fig.subplots_adjust(wspace=0) 
    fig.text(0.5,0.04, r"$\phi$ (deg)", ha='center', va='center') 
    fig.text(0.05,0.5, r"$\psi$ (deg)", ha='center', va='center',rotation='vertical') 
    
    index=0
    for molec in molecList : 
        print index, index/figCols, index%figRows
        ax = axarr[index/figCols,index%figCols]
    
        try : 
            datafile = glob.glob('%s/%s/chi/ramaPhiPsiCNF*.xvg'%(state,molec))[0]
        except IndexError : 
            print "ERROR: Failed to find file '%s/%s/chi/ramaPhiPsiCNF*.xvg"%(state,molec)
            continue 
        resIndex = os.path.splitext(basename(datafile))[0]
        resIndex = resIndex.split("CNF")[1]
    
        datafile2 = glob.glob("%s/GFP_WT/chi/ramaPhiPsi???%s.xvg"%(state,resIndex))[0]
    
        try : 
            data = np.genfromtxt(datafile,skip_header=41) 
        except IOError : 
            print "%s not found!"%datafile
    
        try : 
            data2 = np.genfromtxt(datafile2,skip_header=41) 
        except IOError : 
            print "%s not found!"%datafile2
    
        ax.scatter(data[:,0],data[:,1],s=0.1,color='b') 
        ax.scatter(data2[:,0],data2[:,1],s=0.1,color='g') 
    
        ax.set_title(molec) 
        ax.set_xlim(-180,180)
        ax.set_ylim(-180,180) 
    
        index +=1
    fig.savefig('figures/%s_ramaPhiPsi.png'%(state)[0],format='png') 
    
    
    fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
    fig.subplots_adjust(wspace=0) 
    fig.text(0.5,0.04, r"$\chi_1$ (deg)", ha='center', va='center') 
    fig.text(0.05,0.5, r"$\chi_2$ (deg)", ha='center', va='center',rotation='vertical') 
    
    index=0
    for molec in molecList : 
        print index, index/figCols, index%figRows
        ax = axarr[index/figCols,index%figCols]

        try : 
            datafile = glob.glob('%s/%s/chi/ramaPhiPsiCNF*.xvg'%(state,molec))[0]
        except IndexError : 
            print "ERROR: Failed to find file '%s/%s/chi/ramaPhiPsiCNF*.xvg"%(state,molec)
            continue 
        resIndex = os.path.splitext(basename(datafile))[0]
        resIndex = resIndex.split("CNF")[1]
    
        try : 
            data = np.genfromtxt(datafile,skip_header=41) 
        except IOError : 
            print "%s not found!"%datafile
    
        try : 
            datafile2 = glob.glob("%s/GFP_WT/chi/ramaX1X2???%s.xvg"%(state,resIndex))[0]
            data2 = np.genfromtxt(datafile2,skip_header=41) 
        except (IndexError,IOError) : 
            print "%s not found!"%datafile2
    
        ax.scatter(data[:,0],data[:,1],s=0.1,color='b') 
        ax.scatter(data2[:,0],data2[:,1],s=0.1,color='g') 
    
        ax.set_title(molec) 
        ax.set_xlim(-180,180)
        ax.set_ylim(-180,180) 
    
        index +=1
    fig.savefig('figures/%s_ramaX1X2.png'%(state)[0],format='png') 

