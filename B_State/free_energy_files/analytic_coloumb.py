#!/usr/bin/env python 

import numpy as np 
import sys 
import os.path

try : 
    fileName = sys.argv[1] 
    resName  = sys.argv[2] 
    bVstart  = sys.argv[3] 
    bVend    = sys.argv[4] 
except : 
    print "Usage: %s < pqr file > < residue name > <atom1 of bond vector > < atom2 of bond vector> "%sys.argv[0] 
    print "\t   Potential will be calculated at the midpoint of the bond vector." 
    exit() 
if not os.path.isfile(fileName) : 
    print "ERROR: %s not found"%fileName 
    exit() 

with open(fileName) as f : 
    fileData = f.readlines() 

lineData = [] 
for line in fileData :  
    lineData.append(line.split() ) 

atoms = [] 
for line in lineData  :
    if line[0] == 'ATOM' : 
        atoms.append([line[3], line[2], line[8], line[5], line[6], line[7]]) ## Residue, name, Charge, x, y, z 

for index,atom in enumerate(atoms) : 
    if atom[0] == 'CNF' : 
        if atom[1] == bVstart : 
            atom1 = index 
            atom[2] = 0.00    ##Zero charge on atoms in bond vector (otherwise they dominate the RF) 
        if atom[1] == bVend : 
            atom2 = index 
            atom[2] = 0.00    ##Zero charge on atoms in bond vector (otherwise they dominate the RF) 

try : 
    atom1, atom2
except NameError : 
    print "ERROR: Unable to find %s and %s in residue %s"%(bVstart, bVend, resName) 
    exit() 

atoms = np.array(atoms) 
atoms = np.array(atoms[:,2:],dtype='f')  

dum = np.zeros(3) 
F = 0 

k = 1 / ( 4 * 3.14 * 8.85 * 10 ** -12) ## k = 1/(4*pi*Eps0)  N m^2 C^-2
k = k * (1.602 * 10 **-19 )**2         ## k : N m^2 / e^2
k = k * 10**10 * 10**10                ## k : N ang^2 / e^2 
k = k / (4.11 * 10**-11 )              ## k : kbT / A / Ang^2 / e^2

bV = (atoms[atom1][1:] - atoms[atom2][1:]) / 2
bMP = bV + atoms[atom2][1:]
#print atoms[atom1] 
#print atoms[atom2]
#print bMP

bL = np.sqrt(np.dot(bV, bV) ) 

for atom in atoms : 
    coord = atom[1:] 
    q = atom[0] 
    r = bMP - coord 
    magR = np.sqrt(r[0]**2 + r[1]**2 +r[2]**2)

    F += k*q*r/ magR**3

#print F   ## Potential for charge e, is identical in units of kbT/eA 
print np.dot(F, bV) 
#print np.sqrt(F[0]**2 + F[1]**2 + F[2]**2 ) 

