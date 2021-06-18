#!/usr/bin/env python
import sys
import math
import re
from lxml import etree
from numpy import *
from time import time
import numpy as np
import csv
import os

# Checks whether is running in windows or in other OS
is_windows = sys.platform.startswith('win')
if is_windows:
    slash="\\"
else:
    slash="/"



#Defines the projectname with and the full path to it (the path are the directories where input and output files will reside
#Default values:
inpath = "."
outpath = "."
importfile = ""
projectname = ""
sgbsdenname = ""
pfilename = ""
nbasis=0

#Reads the arguments passed to this interface
argc = len(sys.argv)
#sys.stdout.write("argc = "+str(argc)+"\n")
if sys.argv[1] == '-h':
    sys.stdout.write( "Input should be a MOPAC auxfile of the form <filename>.aux \n")
    quit()
if argc < 2:
    importfile = raw_input('Filename (with or without the extension  .aux): ')
else:
    importfile = sys.argv[1]
    
if (len(importfile) < 4) | (importfile[-4:] != ".aux"):
    projectname = importfile
    importfile = importfile + ".aux"
else:
    projectname = importfile[0:-4]
    
if argc > 2:
    inpath = sys.argv[2]
    outpath = inpath
    if argc > 3:
        outpath = sys.argv[3]
        if argc > 4:
            projectname = sys.argv[4]
if (inpath[-1] != slash):
    inpath += slash
if (outpath[-1] != slash):
    outpath += slash
outputfile = outpath+projectname + "-MOPAC_aux_interface.out"
sys.stdout.write("outputfile = "+ str(outputfile) + "\n")
fout=open(outputfile,'w')
#fout.write("argc = "+ str(argc)+"\n")
#fout.write("argv = ")
#for i in sys.argv:
#    fout.write(i+", ")
#fout.write("\n")
#fout.write("outputfile = "+ outputfile+"\n")
#auxfilename = sys.argv[1]                        # Input Auxfile
#if (len(auxfilename) < 4) | (auxfilename[-4:] != ".aux"):
    #auxfilename = auxfilename + ".aux"
    
auxfilename = inpath + importfile
#sys.stdout.write("auxfilename = "+ auxfilename + "\n")
#fout.write("auxfilename = "+ auxfilename+"\n")
auxfile = open(auxfilename,'r')
sgbsdenname = outpath+projectname+".sgbsden"
#sys.stdout.write( "sgbsdenname = " + sgbsdenname + "\n")
#fout.write("sgbsdenname = "+ sgbsdenname+"\n")
sgbsdenfile = open(sgbsdenname,'w')
pfilename = outpath+projectname+".lowerP"
#fout.write("pfilename = "+ pfilename+"\n")
pfile   = open(pfilename,'w+')            # Density matrix lower triangle in sgsden form



PERIODICTABLE = [' H', 'He', 'Li', 'Be', ' B', ' C', ' N', ' O', ' F', 'Ne',
        'Na', 'Mg', 'Al', 'Si', ' P', ' S', 'Cl', 'Ar', ' K', 'Ca', 'Sc', 'Ti'
        , ' V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 
        'Se', 'Br', 'Kr', 'Rb', 'Sr', ' Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh'
        , 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', ' I', 'Xe', 'Cs', 'Ba', 
        'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er'
        , 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', ' W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 
        'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa'
        , ' U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Mi', 'XX', '+3', '-3', 'Cb',
        '++', ' +', '--', ' -', 'Tv']

def atomic_number(element):
    anum = 1
    for E in PERIODICTABLE:
        if E.strip() == element:
            return anum
        else:
            anum += 1

def type2lm(orbtype):
    if orbtype == 'S': return  (0,0)
    if orbtype == 'PX': return (0)
    if orbtype == 'PY': return (0)
    if orbtype == 'PZ': return (0)
    
def IndicesToLower(u,v): 
    i = u+1
    j = v+1 
    if i >= j:
        return int((i*(i-1)/2.) + j)-1
    else:
        return int((j*(j-1)/2.) + i)-1
        
def LowerToIndices(w):
    k = w+1
    i = int(floor((1+sqrt(8.*k-7))/2.))
    j = int(floor(k - (i*(i-1))/2.))
    return (i-1,j-1)
        
#  THE MAIN PROGRAM
#  ================

fout.write("MOPAC_interface output\n")
fout.write("==========================\n\n")
fout.write("inpath = "+inpath+"\n")
fout.write("outpath = "+outpath+"\n")
fout.write("PROJECT NAME: "+projectname+"\n")
fout.write("auxfilename = "+auxfilename+"\n")
fout.write("sgbsdenfile = "+sgbsdenname+"\n")
fout.write("pfile = "+pfilename+"(removed) \n")


INITIAL_TIME = time()
#sys.stdout.write( "empieza \n" )
DENSITYLEN = 0
SPH_OFFSET = {'S':0,'PX':+2,'PY':-1,'PZ':-1,'X2':4,'XZ':2,'Z2':0,'YZ':-2,'XY':-4}
#sys.stdout.write( "sys.argv[1] = " + str( sys.argv[1]) + "\n")
#auxfile = file(sys.argv[1],'r')                        # Input Auxfile
#filename = sys.argv[1].split(".")[0]
#sys.stdout.write( "filename = " + filename + "\n")
#sgbsdenfile = file(filename+".sgbsden",'w')
#pfile   = file(filename+".lowerP",'w+')            # Density matrix lower triangle in sgsden form

#sys.stdout.write( "sgbsdenfile = " + sgbsdenfile + "\n")

# FLAGS TO TURN RECORDING ON AND OFF AND TO RECORD THE RIGHT VALUES IN THE RIGHT ARRAYS
getAtomEl = False
getAtomCore = False
getGeometry = False
getOrbindices = False
getOrbtypes = False
getZeta = False
getPQN = False
getDensityTriangle = False

# ARRAYS TO STORE INFORMATION RELEVANT FOR OUTPUT
GEOMETRY = []
ORBINDICES = []
ORBTYPES = []
ZETAE = []
PQN = []
ATOM_CORE = []
ATOM_EL = []
NUMBER_OF_CENTERS = 0 # I could get this by taking the length of the GEOMETRY vector. 

mindex  = 0
P = [None]
knt = 0
for line in auxfile:
    # recording flags management
    if "ATOM_EL" in line: # Atomic number
        getAtomEl = True
        continue
    elif "ATOM_CORE" in line: # Atom Core Formal Charge
        getAtomEl = False
        getAtomCore = True
        continue    
    elif "ATOM_X:ANGSTROMS" in line: # Geometry
        getAtomCore = False
        continue
    elif "AO_ATOMINDEX" in line: # Indices of atoms with respect to atomic orbitals
        Number_of_AOs = int(line.split('[')[1].split(']')[0])
#        getGeometry = False
        getAtomCore = False
        getOrbindices = True
        continue
    elif "ATOM_SYMTYPE" in line: 
        getOrbindices = False
        getOrbtypes = True
        continue
    elif "AO_ZETA" in line: # Zeta exponents
        getOrbtypes = False
        getZeta = True
        continue
    elif "ATOM_PQN" in line: # Principal Quantum Numbers
        getZeta = False
        getPQN = True
        continue
    elif "NUM_ELECTRON" in line: 
        getPQN = False
        continue
    elif "ATOM_X_OPT:ANGSTROMS" in line: # Geometry
        NUMBER_OF_CENTERS = int(line.split('[')[1].split(']')[0]) / 3
#        sys.stdout.write( "NUMBER_OF_CENTERS = " + NUMBER_OF_CENTERS + "\n")
        getGeometry = True
        continue
    elif "ATOM_CHARGES" in line:
        getGeometry = False
    elif "DENSITY_MATRIX" in line: # Lower Triangle of Density Matrix in Non-Canonical Order
        getDensityTriangle = True
        DENSITYLEN = int(line.split('[')[1].split(']')[0])
        continue
    elif "]=" in line:
        getDensityTriangle = False
        continue
    else:
        pass 

    # actual stuff management
    
    if "#" in line:
        pass
    elif getAtomEl:
        elements = line.strip().split()
        for element in elements:
            ATOM_EL.append(atomic_number(element))
#        sys.stdout.write( "ATOM_EL = " + ATOM_EL + "\n")
    elif getAtomCore:
        elements = line.strip().split()
        for element in elements:
            ATOM_CORE.append(int(element))
    elif getOrbindices:
        elements = line.strip().split()
        for element in elements:
            ORBINDICES.append(int(element))
    elif getOrbtypes:
        elements = line.strip().split()
        for element in elements:
            ORBTYPES.append(element)
    elif getZeta:
        elements = line.strip().split()
        for element in elements:
            ZETAE.append(float(element))
    elif getPQN:
        elements = line.strip().split()
        for element in elements:
            PQN.append(int(element))
    elif getGeometry:
        point    = list(map(lambda x: float(x), line.strip().split()))
        GEOMETRY.append(point)
#        sys.stdout.write( "GEOMETRY = " + str( GEOMETRY ) + "\n")
    elif getDensityTriangle:
        elements = line.strip().split()
        for element in elements:
            old_i,old_j = LowerToIndices(knt)
            new_i = old_i+SPH_OFFSET[ORBTYPES[old_i]]
            new_j = old_j+SPH_OFFSET[ORBTYPES[old_j]]
            new_k = IndicesToLower(new_i,new_j)
            if len(P) == 1: mindex = knt
            current_index = new_k - mindex
            while(len(P)-1) < current_index: P.append(None)
            line = "%22.15e " % float(element)
            if ((new_k+1)%8==0): line += "\n"
            P[current_index] = line
            if None not in P: 
                for element in P: pfile.write(element)
                P = [None]
                mindex = 0
            knt += 1
    else: 
        pass

del P
pfile.seek(0)

if 1 == 2: # NEVER RECORD THESE THINGS, GOTTA FIND A WAY TO TURN THESE ON AND OFF BASED ON COMMAND-LINE
    # Defining the energy vectors
    EIGENVALUES = []
    i = auxfile.index(filter(lambda x: (("EIGENVALUES" in x) or ("LMO_ENERGY" in x)), auxfile)[0])
    j = auxfile.index(filter(lambda x: "MOLECULAR_ORBITAL_OCCUPANCIES" in x,auxfile)[0])
    eigenvals = auxfile[i+1:j]
    for line in eigenvals:
        elements = line.strip().split()
        for element in elements:
            EIGENVALUES.append(float(element))

    # Defining occupancies vectors
    OCCUPANCIES = []
    i = auxfile.index(filter(lambda x: "MOLECULAR_ORBITAL_OCCUPANCIES" in x, auxfile)[0])
    j = auxfile.index(filter(lambda x: "CPU_TIME:SECONDS" in x,auxfile)[0])
    occupancies = auxfile[i+1:j]
    for line in occupancies:
        elements = line.strip().split()
        for element in elements:
            OCCUPANCIES.append(float(element))


    # Defining coefficient matrices
    C  = []
    MO = []
    i  = auxfile.index(filter(lambda x: (("EIGENVECTORS" in x) or ("LMO_VECTORS" in x)), auxfile)[0])
    j  = auxfile.index(filter(lambda x: "DENSITY_MATRIX" in x,auxfile)[0])
    coefficients = auxfile[i+1:j]
    for line in coefficients:
        elements = line.strip().split()
        for element in elements:
            if len(MO) >= Number_of_AOs:
                C.append(MO)
                MO = []
            else:
                MO.append(float(element))

    # Writting file with occupancy, coefficient and molecular orbital matrices. 
    i = 0
    for MO in C:
        line =  "%22.15e  %22.15e" % (OCCUPANCIES[i],EIGENVALUES[i])
        cfile.write(line)
        for k in MO:
            line = " %22.15e" % k
            cfile.write(line)    
        i += 1
        cfile.write('\n')


# Compiling all the information about orbitals in a single master array of strings
GEOMETRYLINES = []
SHELLLINES = []
ORBS = {'S':0,'P':1,'X':2,'Z':2,'Y':2} # A way to link the strings of atom types with their "L" quantum numbers
ORBINDICES+=[0] # It needs this last pseudo-atom-index in order for the loop not to overflow the array
orb_i = 0
#sys.stdout.write( "NUMBER_OF_CENTERS = " + str( NUMBER_OF_CENTERS ) + "\n")
#sys.stdout.write( "GEOMETRY = " + str( GEOMETRY ) + "\n")
chargescenter = [0,0,0]
totalcharge = 0.
for i in range(int(NUMBER_OF_CENTERS)):
    for j in range(3):
        chargescenter[j] = chargescenter[j] + GEOMETRY[i][j] * ATOM_EL[i]
    totalcharge += ATOM_EL[i]
for i in range(3):
    chargescenter[i] = chargescenter[i] / totalcharge
#sys.stdout.write( "totalcharge = " + str( totalcharge ) + "\n")
#sys.stdout.write( "chargescenter = " + str( chargescenter ) + "\n")

for i in range(int(NUMBER_OF_CENTERS)):
    center = GEOMETRY[i]
    x      = (center[0]-chargescenter[0]) * 1.889725989 # conversion to Bohr
    y      = (center[1]-chargescenter[1]) * 1.889725989
    z      = (center[2]-chargescenter[2]) * 1.889725989
    atom_core = ATOM_CORE[i]
    atom_el = ATOM_EL[i]
    zeta = 0.0
    nshells = 0
    n = PQN[orb_i]
    while (ORBINDICES[orb_i]-1 == i):        # Cycle through all the orbitals that have the same center. Orbitals and Atoms are (...)    
                                            # (...) ordered such that orb_i never needs to retrocede.
        if ZETAE[orb_i] == zeta:            # If the Zetas are equal, this shell has been accounted for already
            pass
        else:
            ll = ORBS[ORBTYPES[orb_i][0]]
            zeta = ZETAE[orb_i]
            shelline = "%2i %2i %22.15e \n" % (n,ll,zeta)
            SHELLLINES.append(shelline)
            nshells += 1
        orb_i += 1
    geoline = "%22.15e %22.15e %22.15e %2i %2i \n" % (x,y,z,atom_el,nshells)
#    geoline = "%22.15e %22.15e %22.15e %3i %2i %2i \n" % (x,y,z,atom_el,nshells,atom_core)
    GEOMETRYLINES.append(geoline)    
        

# Cleaning memory
del PQN
del ORBINDICES
del GEOMETRY
OUTLINES = [str(NUMBER_OF_CENTERS)+" \n"] + GEOMETRYLINES + SHELLLINES
del GEOMETRYLINES
del SHELLLINES
sgbsdenfile.writelines(OUTLINES)
for line in pfile: sgbsdenfile.write(line)
pfile.close()
sgbsdenfile.close()
os.remove(pfilename)

fout.write("Timing: "+str(time() - INITIAL_TIME)+" seconds. \n")
sys.stdout.write( "It all took: "+str(time() - INITIAL_TIME)+" seconds. \n" )
quit()
