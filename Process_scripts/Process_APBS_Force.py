import glob 
import numpy as np 
import matplotlib.pyplot as plt
from sys import exit
import os 
from scipy.stats import linregress

from matplotlib import rc_file 

save_dir='figures/force_calc'
experimental_data='Exp_data/sasa_abs_data.dat'

figCols=3
figRows=3

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

nameToColorKeys = {
        "GFP_F145X":'#A93226',
        "GFP_M218X":'r',
        "GFP_F165X":'#E67E22',
        "GFP_F114X":'y',
        "GFP_N198X":'#2ECC71',
        "GFP_Y143X":'g',
        "GFP_N212X":'c',
        "GFP_D117X":'b',
        "GFP_Y92X":'m'
        }

fig, axarr = plt.subplots(1,1) 
fig.set_size_inches(3.5,3) 
fig.subplots_adjust(hspace=0.5,left=0.2,top=0.9,bottom=0.1) 

ax = axarr

force_files = glob.glob('B_State/*/APBS_force_calc/rxn_field.out') 

avgAccum = [] 
expAccum = [] 
for fFile in force_files : 
    molec = fFile.split('/')[1]
    print molec 

    try : 
        rxn_field = np.genfromtxt(fFile) 
        coloumb_field = np.genfromtxt('B_State/%s/APBS_force_calc/coloumb_field.out'%molec) 
    except : 
        print "%s file failed to import. Skipping."%fFile 
        continue 

    data = rxn_field + coloumb_field
    #data = coloumb_field
    #data = rxn_field
    print len(data),
    data = data[len(data)/5:]
    print len(data) 

    avg = np.average(data) 
    std = np.std(data) 
    exp = float(mutToExp[molec]) 

    avgAccum.append(avg) 
    expAccum.append(exp) 
    print "%12s%10.1f +/- %3.1f\t%.1f"%(molec, avg, std, exp) 

    ax.errorbar(exp, avg, std, ecolor='k',fmt='none') 
    ax.scatter(exp,avg,marker='D',edgecolor='k',color=nameToColorKeys[molec],zorder=3) 

ax.axvline(2227.5, color='k', linestyle='--') 
ax.axvline(2235.9, color='b', linestyle='--') 
slope, intercept, r_value, p_value, std_error = linregress(expAccum, avgAccum) 
print "*** r = %0.5f***"%(r_value) 
xs = np.linspace(min(expAccum),max(expAccum),100) 
ys = slope * xs + intercept 
ax.plot(xs,ys,label="r = %0.2f ; slope = %0.2f "%(r_value,slope) ) 

slope=1/1.7219 ; intercept = intercept - 1753
ys = slope * xs + intercept 
ax.plot(xs,ys,label=r"slope = %0.2f $(k_b T/e^{-} \AA)/cm^{-1}$"%(slope) ) 
ax.legend(loc=2,fontsize='xx-small') 
    
fig.savefig('%s/APBS_forces.pdf'%save_dir,format='pdf') 
plt.close() 

fig, ax = plt.subplots(1,1) 
for i in range(len(avgAccum)) : 
    molec = force_files[i].split('/')[1]
    avgAccum[i] /= slope
    ax.scatter(expAccum[i],avgAccum[i],color=nameToColorKeys[molec] ) 
    print "%10s  %5.1f  %5.3f"%(molec, expAccum[i],avgAccum[i]) 
intercept = -2235
ax.plot(xs, xs + intercept) 
plt.show() 


fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row')
fig.subplots_adjust(wspace=0)
fig.text(0.5,0.04, r"R$_g$ ($\AA$)", ha='center', va='center')
fig.text(0.03,0.5, r"RMSD ($\AA$)", ha='center', va='center',rotation='vertical')

figD, axarrD = plt.subplots(figRows,figCols,sharex='col',sharey='row')
figD.subplots_adjust(wspace=0)
figD.text(0.5,0.04, r"R$_g$ ($\AA$)", ha='center', va='center')
figD.text(0.03,0.5, r"RMSD ($\AA$)", ha='center', va='center',rotation='vertical')

molList = glob.glob('B_State/*/force_calc') 
for index, mol in enumerate(molList) :
    molec = mol.split('/')[1]
    print index, index/figCols, index%figRows
    ax = axarr[index/figCols,index%figCols]
    axD = axarrD[index/figCols,index%figCols]

    try : 
        rxn_field = np.genfromtxt(fFile) 
        coloumb_field = np.genfromtxt('B_State/%s/APBS_force_calc/coloumb_field.out'%molec) 
    except : 
        print "%s file failed to import. Skipping."%fFile 
        continue 

    data = rxn_field + coloumb_field
    data = data[len(data)/5:]

    avg, firstD = [],[]
    for i in range(0,len(data[:]) ) :
        avg.append(np.average(data[:i]))
    for i in range(len(avg) - 1) :
        firstD.append(avg[i+1] - avg[i])
    firstD = np.array(firstD)

    ax.plot(avg,color='b')
    ax.set_title(molec.split('_')[1]) 
    axD.plot(firstD,color='b') 
    
    axD.set_ylim([-1,1]) 

fig.savefig('%s/APBS_forces-time.png'%save_dir,format='png') 
figD.savefig('%s/APBS_forces-timeD.png'%save_dir,format='png') 
