import glob 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm
import os
from sys import exit
from scipy.optimize import curve_fit

figCols=3
figRows=3

if not os.path.isdir('figures') : 
    os.mkdir('figures') 


fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0,hspace=0.30) 
fig.subplots_adjust(left=0.1,right=0.8,bottom=0.1,top=0.9) 
fig.text(0.5,0.04, "Time (ns)", ha='center', va='center') 

def double_exp(x, a, b, c, d) : 
    return a * np.exp(b * x) + c * np.exp(d * x) 
def single_exp(x, a, b, c, d) : 
    return a * np.exp(b * x) 
def linear_fit(x, m, b) : 
    return m*x+b

molecProtRates = {} 
molecWatRates= {} 
xmax = 0 
for donor in ["water"] :#, "protein"] : 
    for temp in [290, 300, 310] :     
        if donor == "water" : 
            col = 4
            if temp == 290 : 
                datafiles = glob.glob('B_State/*/hbond_with_ca_290/persistent.xvg') 
                color = 'navy' 
            elif temp == 310 : 
                datafiles = glob.glob('B_State/*/hbond_with_ca_310/persistent.xvg') 
                color = 'blue'
            else : 
                datafiles = glob.glob('B_State/*/hbond_with_ca/persistent.xvg')    
                color = 'deepskyblue'
        #elif donor == "protein" : 
        #    col = 6 
        #    if temp == 290 : 
        #        datafiles = glob.glob('B_State/*/hbond_with_ca_290/persistent.xvg') 
        #        color = 'darkgreen' 
        #    elif temp == 310 : 
        #        datafiles = glob.glob('B_State/*/hbond_with_ca_310/persistent.xvg') 
        #        color = 'limegreen'
        #    else : 
        #        datafiles = glob.glob('B_State/*/hbond_with_ca/persistent.xvg')    
        #        color = 'green'
        for index,file in enumerate(datafiles) : 
            try : 
                data = np.genfromtxt(file,skip_header=24) 
                x,y = data[:,0], data[:,col]
                x /= 1000
                x,y = x[1:], y[1:] 
                y = np.trim_zeros(y)  ##Trim out times where zero hbonds lasted
                x = x[:len(y)]          ## otherwise exponential fit is difficult 
                y = np.log(y) 
            except : 
                print "Import failed for %s" %file
                continue 

            #if (np.max(data[:,0]) > xmax) : xmax = np.max(data[:,0]) 

            ax = axarr[index/figCols,index%figCols]

            ax.scatter(x,y, marker='o', s=1.0, c= color, edgecolors='none' ) 
            #ax.set_xscale('log') 
            #ax.set_yscale('log') 
            ax.set_title(file.split('/')[1].split('_')[1]) 
            #ax.set_ylim([0.8,1000]) 

            #try : 
            #    popt, pcov = curve_fit(double_exp, x, y, p0=(1, -1, 1, -1))
            #    p1, p2, p3, p4 = popt[:]
            #    #print popt
            #    ax.plot(x, double_exp(x, p1, p2, p3, p4), c = color, linewidth=1, label=donor+"-"+str(temp)+" K") 
            #except : 
            #    try : 
            #        popt, pcov = curve_fit(double_exp, x, y, p0=(1, -1, 1, -1))
            #        p1, p2, p3, p4 = popt[:]
            #        #print popt
            #        ax.plot(x, double_exp(x, p1, p2, p3, p4), c = color) 
            #    except : 
            #        continue 
            #molec = file.split('/')[1].split('_')[1]
            #if popt[0] > popt[2] : rate = popt[1] 
            #else : rate = popt[3]
            #
            #if donor == "protein" : 
            #    if not molec in molecProtRates : 
            #        molecProtRates[molec] = [[temp,rate]]
            #    else : 
            #        molecProtRates[molec].append([temp,rate]) 
            #elif donor == "water" : 
            #    if not molec in molecWatRates : 
            #        molecWatRates[molec] = [[temp,rate]]
            #    else : 
            #        molecWatRates[molec].append([temp,rate]) 

            #ax.loglog(data[:,0],data[:,4],color='c',label='Nearby water') 
            #ax.loglog(data[:,0],data[:,5],color='b', label='HBonding water') 
            #ax.loglog(data[:,0],data[:,6],color='g',label='Hbonding protein') 


for index,file in enumerate(datafiles) : 
    ax = axarr[index/figCols,index%figCols]
    ax.set_xlim(0.004,0.15) 

fig.legend(loc='center right', fontsize='x-small') 
fig.savefig('figures/Exponential_fits.pdf',format='pdf') 
plt.close() 

#print molecProtRates
#print molecWatRates
for molec in molecWatRates :
    data = molecWatRates[molec]
    data = np.array(data,dtype='float') 
    x = 1/data[:,0]
    y = np.log(-1*data[:,1]) 
    #print y 
    #print data
    #x = 1/data[:,0]
    #y = np.log(data[:,1]) 
    #print data 
    plt.plot(x,y,label=molec) 

plt.legend(loc=2) 
plt.xlabel("1/T(K)") 
plt.ylabel("ln(k)") 
plt.savefig('figures/Arrhenius.pdf',format='pdf') 
