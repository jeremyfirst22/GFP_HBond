import glob 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm
import os
from sys import exit
from scipy.optimize import curve_fit
from matplotlib import rc_file
from scipy.stats import linregress

rcFile = 'rc_files/paper.rc'
rc_file(rcFile) 

figCols=2
figRows=2

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

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
        #"Y92X"]#,
        "F114X",
        "D117X",
        #"Y143X",
        #"F145X"]#,
        #"F165X",
        "N198X",
        "N212X"]#,
        #"M218X"]


fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
#fig, axarr = plt.subplots(1,1)
fig.subplots_adjust(wspace=0,hspace=0.30) 
fig.subplots_adjust(left=0.1,right=0.8,bottom=0.1,top=0.9) 
fig.text(0.5,0.04, "Time (ns)", ha='center', va='center') 
fig.text(0.02,0.50, "Number of hbonds", ha='center', va='center',rotation='vertical') 

def double_exp(x, a, b, c, d) : 
    return a * np.exp(b * x) + c * np.exp(d * x) 
def single_exp(x, a, b) : 
    return a * np.exp(b * x) 

molecProtRates = {} 
molecWatRates= {} 
xmax = 0 
for index, molec in enumerate(molecList) : 
    #for donor in ["water", "protein"] : 
    for donor in ["water"] : 
        for temp in [290, 300, 310] :     
            if donor == "water" : 
                col = 4
                if temp == 290 : 
                    file = 'B_State/GFP_%s/hbond_with_ca_290/persistent.xvg'%molec  
                    color = 'navy' 
                elif temp == 310 : 
                    file = 'B_State/GFP_%s/hbond_with_ca_310/persistent.xvg'%molec  
                    color = 'blue'
                else : 
                    file = 'B_State/GFP_%s/hbond_with_ca/persistent.xvg'%molec  
                    color = 'deepskyblue'
            elif donor == "protein" : 
                col = 6 
                if temp == 290 : 
                    file = 'B_State/GFP_%s/hbond_with_ca_290/persistent.xvg'%molec  
                    color = 'darkgreen' 
                elif temp == 310 : 
                    file = 'B_State/GFP_%s/hbond_with_ca_310/persistent.xvg'%molec  
                    color = 'limegreen'
                else : 
                    file = 'B_State/GFP_%s/hbond_with_ca/persistent.xvg'%molec  
                    color = 'green'

            try : 
                data = np.genfromtxt(file,skip_header=24) 
                x,y = data[:,0], data[:,col]
                x /= 1000
                x,y = x[1:], y[1:] 
                y = np.trim_zeros(y)  ##Trim out times where zero hbonds lasted
                x = x[:len(y)]          ## otherwise exponential fit is difficult 
            except : 
                print "Import failed for %s" %file
                continue 
    
            #if (np.max(data[:,0]) > xmax) : xmax = np.max(data[:,0]) 
    
            ax = axarr[index/figCols,index%figCols]
            #ax = axarr
    
            ax.scatter(x,y, marker='o', s=2.0, c= color, edgecolors='none' ) 
            #ax.set_xscale('log') 
            ax.set_yscale('log') 
            ax.set_title(molec) 
            ax.set_ylim([0.8,1000]) 
    
            if molec == "F165X" or (donor == "protein" and molec == "Y92X") or (donor == "water" and molec == "F145X"): 
                popt, pcov = curve_fit(double_exp, x, y, p0=(1, -1, 1, -1))
                p1, p2, p3, p4 = popt[:]
                print temp, popt
                ax.plot(x, double_exp(x, p1, p2, p3, p4), c = color, linewidth=1, label=donor+"-"+str(temp)+" K") 

                if popt[0] > popt[2] : rate = popt[1] 
                else : rate = popt[3]
            #except : 
            else:
                #try : 
                popt, pcov = curve_fit(single_exp, x, y, p0=(1, -1, ))
                p1, p2, = popt[:]
                #print popt
                ax.plot(x, single_exp(x, p1, p2), c = color, linewidth=1, label=donor+"-"+str(temp)+" K")

                rate = popt[1]
                #except : 
                #    print "Warning: Single exponential failed for %s %s"%(donor, molec)
                #    continue 
            
            if donor == "protein" : 
                if not molec in molecProtRates : 
                    molecProtRates[molec] = [[temp,rate]]
                else : 
                    molecProtRates[molec].append([temp,rate]) 
            elif donor == "water" : 
                if not molec in molecWatRates : 
                    molecWatRates[molec] = [[temp,rate]]
                else : 
                    molecWatRates[molec].append([temp,rate]) 
    
            #ax.loglog(data[:,0],data[:,6],color='g',label='Hbonding protein') 


for index,file in enumerate(molecList) : 
    ax = axarr[index/figCols,index%figCols]
    #ax = axarr
    ax.set_xlim(0.004,0.01 ) 

fig.legend(loc=1) 
fig.savefig('figures/Exponential_fits.pdf',format='pdf') 
plt.close() 

#print molecProtRates
#print molecWatRates
for molec in molecWatRates :
    data = molecWatRates[molec]
    data = np.array(data,dtype='float') 
    x = 1/data[:,0]
    y = np.log(-1*data[:,1]) 

    slope, intercept, r_value, p_value, std_error = linregress(x,y) 
    #print y 
    #print data
    #x = 1/data[:,0]
    #y = np.log(data[:,1]) 
    #print data 
    plt.scatter(x,y,label=molec,color=nameToColorKeys[molec]) 
    plt.plot(x,x*slope+intercept,color=nameToColorKeys[molec]) 

plt.xlim(0.00322,0.00346) 
plt.legend(loc=3) 
plt.xlabel("1/T(K)") 
plt.ylabel("ln(k)") 
plt.savefig('figures/Arrhenius_water.pdf',format='pdf') 
plt.close() 

for molec in molecProtRates :
    data = molecProtRates[molec]
    data = np.array(data,dtype='float') 
    x = 1/data[:,0]
    y = np.log(-1*data[:,1]) 

    slope, intercept, r_value, p_value, std_error = linregress(x,y) 
    #print y 
    #print data
    #x = 1/data[:,0]
    #y = np.log(data[:,1]) 
    #print data 
    plt.scatter(x,y,label=molec,color=nameToColorKeys[molec]) 
    plt.plot(x,x*slope+intercept,color=nameToColorKeys[molec]) 


plt.xlim(0.00322,0.00346) 
plt.legend(loc=2) 
plt.xlabel("1/T(K)") 
plt.ylabel("ln(k)") 
plt.savefig('figures/Arrhenius_protein.pdf',format='pdf') 
