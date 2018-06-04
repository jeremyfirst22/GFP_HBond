import glob 
import numpy as np 
import matplotlib.pyplot as plt
from sys import exit
import os 
from matplotlib import rc_file
from scipy.stats import linregress

from matplotlib import rc_file 

save_dir='figures/force_calc'
experimental_data='Exp_data/sasa_abs_data.dat'

rcFile = 'rc_files/paper.rc'

figCols=4
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

rc_file(rcFile) 

#fig, axarr = plt.subplots(1,1) 
fig, ax = plt.subplots(1,1,figsize=(4.3,3) ) 
fig.subplots_adjust(left=0.15,right=0.88,bottom=0.15)
#fig.set_size_inches(3.5,6) 
#fig.subplots_adjust(hspace=0.5,left=0.2,top=0.9,bottom=0.1) 
#for index, force_type in enumerate(['external_field','solvent_rxn_field','total_field']) : 
for index, force_type in enumerate(['external_field']) : 
    #ax = axarr[index] 
    #ax.set_title(force_type ) 

    force_files = glob.glob('B_State/*/force_calc/*.%s.projected.xvg'%force_type) 

    avgAccum = [] 
    expAccum = [] 
    for fFile in force_files : 
        molec = fFile.split('/')[3].split('.')[0]  
        print molec 

        try : 
            data = np.genfromtxt(fFile) 
        except : 
            print "%s file failed to import. Skipping."%fFile 
            continue 
    
        avg = np.average(data) 
        std = np.std(data) 
        exp = float(mutToExp[molec]) 

        avgAccum.append(avg) 
        expAccum.append(exp) 
        print "%20s%12s%10.1f +/- %3.1f\t%.1f"%(force_type, molec, avg, std, exp) 

        ax.errorbar(exp, avg, std, label=molec.split('_')[1],marker='o',color=nameToColorKeys[molec],markeredgecolor='none',capsize=3,linestyle='none') 
        #ax.scatter(exp,avg,marker='D',color=nameToColorKeys[molec],zorder=3) 


    slope, intercept, r_value, p_value, std_error = linregress(expAccum, avgAccum) 
    print "***%20s: r = %0.5f***"%(force_type, r_value) 
    xs = np.linspace(min(expAccum),max(expAccum),100) 
    ys = slope * xs + intercept 

    first_legend = plt.legend(bbox_to_anchor=(1.38,0.95),numpoints=1,edgecolor='k') 
    plt.gca().add_artist(first_legend)

    bestfit, = plt.plot(xs,ys,label="r = %0.3f"%r_value,color='k',zorder=1)

    fig.text(0.62,0.82,r"r = %0.3f"%r_value)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width *0.85, box.height])

    
plt.xlabel(r"$\tilde{\nu}$ (cm$^{-1}$)") 
plt.ylabel(r"Calculated Field (k$_b$T/$e^- \rm{\AA})$") 

fig.savefig('%s/forces.pdf'%save_dir,format='pdf') 
#
#fig, axarr = plt.subplots(2,1) 
#fig.set_size_inches(3.5,6) 
#fig.subplots_adjust(hspace=0.5,left=0.2,top=0.9,bottom=0.1) 
#for index, force_type in enumerate(['external_field','total_field']) : 
#    ax = axarr[index] 
#    ax.set_title(force_type+'-solvent_rxn' ) 
#
#    force_files = glob.glob('B_State/*/force_calc/*.%s.projected.xvg'%force_type) 
#
#    avgAccum = [] 
#    expAccum = [] 
#    for fFile in force_files : 
#        molec = fFile.split('/')[3].split('.')[0]  
#        print molec 
#
#        try : 
#            data = np.genfromtxt(fFile) 
#        except : 
#            print "%s file failed to import. Skipping."%fFile 
#            continue 
#        try : 
#            data2 = np.genfromtxt('B_State/%s/force_calc/%s.solvent_rxn_field.projected.xvg'%(molec,molec)) 
#        except : 
#            print "Solvent_rxn_field for %s file failed to import. Skipping."%molec
#            continue 
#        
#        assert len(data) == len(data2) 
#
#        data = data - data2
#
#        data = data[2500:] ##Discard first 10ns as equilibration time. 2500 frames*4ps/frame = 10 ns 
#
#        avg = np.average(data) 
#        std = np.std(data) 
#        exp = float(mutToExp[molec]) 
#
#        avgAccum.append(avg) 
#        expAccum.append(exp) 
#        print "%20s%12s%10.1f +/- %3.1f\t%.1f"%(force_type, molec, avg, std, exp) 
#
#        ax.errorbar(exp, avg, std, ecolor='k',fmt='none') 
#        ax.scatter(exp,avg,marker='D',edgecolor='k',color=nameToColorKeys[molec],zorder=3) 
#
#
#    slope, intercept, r_value, p_value, std_error = linregress(expAccum, avgAccum) 
#    print "***%20s: r = %0.5f***"%(force_type, r_value) 
#    xs = np.linspace(min(expAccum),max(expAccum),100) 
#    ys = slope * xs + intercept 
#    ax.plot(xs,ys,label="r = %0.2f"%r_value) 
#    ax.legend(loc=2,fontsize='xx-small') 
#    
#fig.savefig('%s/forces-minus-solvent-rxn-field.pdf'%save_dir,format='pdf') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row')
fig.subplots_adjust(wspace=0)
fig.text(0.5,0.04, r"Time (frames)",ha='center', va='center')
fig.text(0.03,0.5, r"$\frac{\partial <F>}{\partial t}$",ha='center', va='center',rotation='vertical')

figD, axarrD = plt.subplots(figRows,figCols,sharex='col',sharey='row')
figD.subplots_adjust(wspace=0)
figD.text(0.5,0.04, r"Time (frames)",ha='center', va='center')
figD.text(0.03,0.5, r"$\frac{\partial <F>}{\partial t}$",ha='center', va='center',rotation='vertical')

molList = glob.glob('B_State/*/force_calc') 
for index, mol in enumerate(molList) :
    molec = mol.split('/')[1]
    print index, index/figCols, index%figRows
    ax = axarr[index/figCols,index%figCols]
    axD = axarrD[index/figCols,index%figCols]

    field='external_field'
    datafile = 'B_State/'+molec+'/force_calc/'+molec+'.%s.projected.xvg'%field
    try :
        data1 = np.genfromtxt(datafile) 
        data1 = data1[2500:]
    except IOError :
        print "No file found for %s"%(datafile)
        continue 

    avg, firstD = [],[]
    for i in range(0,len(data1[:]) ) :
        avg.append(np.average(data1[:i]))
    for i in range(len(avg) - 1) :
        firstD.append(avg[i+1] - avg[i])
    firstD = np.array(firstD)

    ax.plot(avg,color='b')
    axD.plot(firstD,color='b',linewidth=0.1)
    
    axD.set_ylim([-.1,.1]) 

fig.savefig('%s/forces-time.png'%save_dir,format='png') 
figD.savefig('%s/forces-timeD.png'%save_dir,format='png') 
