#  Copyright 2013-2021, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,
#  Guillermo Ramirez, Anmol Kumar, Sachin D. Yeole, Shridhar R. Gadre
# 
#  This file is part of DAM320.
# 
#  DAM320 is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  DAM320 is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with DAM320.  If not, see <http://www.gnu.org/licenses/>.
#
#------------------------------------------------------------------------
#
#    Program for generating grids for plotting orbitals
#
# Version of September 2017
#
#!/usr/bin/env python
import sys
import math
import re
from lxml import etree
from numpy import *
import numpy as np
import csv
import os.path

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
nbasis=0

try: input = raw_input
except NameError: pass

#Reads the arguments passed to this interface
argc = len(sys.argv)
sys.stdout.write( "argc = " + str(argc) + "\n")
if argc < 2:
    importfile = input('Filename (with or without the extension  .xml): ')
else:
    importfile = sys.argv[1]

if (len(importfile) < 4) | (importfile[-4:] != ".xml"):
    projectname = importfile
    importfile = importfile + ".xml"
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
outputfile = outpath+projectname + "-MOLPRO_xml_interface.out"

fout=open(outputfile,'w')

if not (os.path.isfile(inpath+importfile)):
    fout.write(inpath+importfile+" does not exist\n")
    raise Exception(inpath+importfile+" does not exist\n")

#sys.stdout.write( "inpath+importfile = " + inpath+importfile + "\n")
tree = etree.parse(inpath+importfile)
#sys.stdout.write( "tree = " + tree + "\n")
lvec = []
document = tree.getroot()

namespaces={
    'molpro-output': 'http://www.molpro.net/schema/molpro-output',
    'xsd': 'http://www.w3.org/1999/XMLSchema',
    'cml': 'http://www.xml-cml.org/schema',
    'stm': 'http://www.xml-cml.org/schema',
    'xhtml': 'http://www.w3.org/1999/xhtml',
    'xlink': 'http://www.w3.org/1999/xlink'}
    
atomicnumber = {"H": 1, "HE": 2, "LI": 3, "BE": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "NE": 10, "NA": 11, "MG": 12,
"AL": 13, "SI": 14, "P": 15,"S": 16, "CL": 17, "AR": 18, "K": 19, "CA": 20, "SC": 21, "TI": 22, "V": 23, "CR": 24, "MN": 25,
"FE": 26, "CO": 27, "NI": 28, "CU": 29,"ZN": 30, "GA": 31, "GE": 32, "AS": 33, "SE": 34, "BR": 35, "KR": 36, "RB": 37,
"SR": 38, "Y": 39, "ZR": 40, "NB": 41, "MO": 42, "TC": 43, "RU": 44, "RH": 45, "PD": 46, "AG": 47, "CD": 48, "IN": 49, 
"SN": 50, "SB": 51, "TE": 52, "I": 53, "XE": 54, "CS": 55, "BA": 56, "LA": 57, "CE": 58, "PR": 59, "ND": 60, "PM": 61,
"SM": 62, "EU": 63, "GD": 64, "TB": 65, "DY": 66, "HO": 67, "ER": 68, "TM": 69, "YB": 70, "LU": 71, "HF": 72, "TA": 73,
"W": 74, "RE": 75, "OS": 76, "IR": 77, "PT": 78, "AU": 79, "HG": 80, "TL": 81, "PB": 82, "BI": 83, "PO": 84, "AT": 85,
"RN": 86, "FR": 87, "RA": 88, "AC":89, "TH": 90, "PA": 91, "U": 92, "NP": 93, "PU": 94, "AM": 95, "CM": 96, "BK": 97,
"CF": 98, "ES": 99, "FM": 100, "MD": 101, "NO": 102, "LR": 103}


def writeBasis(molecule): # Prints basis set
    sys.stdout.write( "molecule = " + str(molecule) + "\n")
    atoms = molecule.xpath('cml:atomArray/cml:atom',namespaces=namespaces)    # Kept for compatibility with old molpro xml files
    if len(atoms) == 0:
        atoms = molecule.xpath('cml:molecule/cml:atomArray/cml:atom',namespaces=namespaces)
    if len(atoms) == 0:
        fout.write("Error reading atoms from cml:molecule/cml:atomArray/cml:atom\n")
        raise Exception('Error reading atoms from cml:molecule/cml:atomArray/cml:atom')
    sys.stdout.write( "molecule.xpath = " + str(molecule.xpath) + "\n")
    sys.stdout.write( "atoms = " + str(atoms) + "\n")
    basisSets = molecule.xpath('molpro-output:basisSet[@id="ORBITAL"]',namespaces=namespaces)
    sys.stdout.write( "basisSets = " + str(basisSets) + "\n")
    if len(basisSets) != 1:
        raise Exception('something is wrong: there should be just one orbital basisSet')
    basisSet = basisSets[0]
    #sys.stdout.write(('Basis length ='+str(basisSet.get('length'))+' angular ='+str(basisSet.get('angular'))+' type =',str(basisSet.get('type'))+' groups ='+str(basisSet.get('groups')) + "\n"))
    basisGroups = basisSet.xpath('molpro-output:basisGroup',namespaces=namespaces)
    # now walk through basis set
    iatom=0
    datosbase = []
    nbas = 0
    for atom in atoms:
        kntcontr = 0
        baseatom = []
        query = 'molpro-output:association/molpro-output:atoms[@xlink:href[contains(.,"@id=\'' + atom.get('id') + '\'")]]/..'
        basisGroupAssociation = basisSet.xpath(query,namespaces=namespaces)
        if len(basisGroupAssociation) != 1:
            raise Exception('something is wrong: there should be a unique association of an atom with a basis set')
        bases=basisGroupAssociation[0].xpath('molpro-output:bases',namespaces=namespaces)
        if len(bases) != 1:
            raise Exception('something is wrong: there should be a bases node in association')
        basesString=bases[0].get('{'+namespaces['xlink']+'}href')
        basesString=basesString[basesString.find('basisGroup['):]
        basesString=basesString[basesString.find('[')+1:].lstrip()
        basesString=basesString[:basesString.find(']')].rstrip()
        basesString=basesString.replace('or','').replace('\n','').replace("'",'')
        list = basesString.split("@id=");
        #sys.stdout.write( "list = " + list + "\n")
        for item in list:
            item=item.lstrip().rstrip()
            if item.isalnum() :
                basisGroup=basisSet.xpath('molpro-output:basisGroup[@id="'+item+'"]',namespaces=namespaces)[0]
                lquant =int(basisGroup.get('minL'))
                if lquant > 6:
                    raise Exception("Sorry, I was too lazy to write this for j basis functions and higher")
                if lquant != int(basisGroup.get('maxL')):
                    raise Exception("This program cannot handle multiple-angular momentum sets")
                #sys.stdout.write(str(basisGroup.get('id'))+' '+str(basisGroup.get('minL'))+' '+str(basisGroup.get('maxL'))+' '+str(basisGroup.get('primitives'))+' '+str(basisGroup.get('angular'))+' '+str(basisGroup.get('contractions')= + "\n"))
                alpha=np.float64(re.sub(' +',' ',basisGroup.xpath('molpro-output:basisExponents',namespaces=namespaces)[0].text.replace('\n','').lstrip().rstrip()).split(" "))
                #sys.stdout.write(alpha + "\n"))   # Prints the exponents
                #first evaluate all the primitives for this shell on the grid
                lcontr = basisGroup.get('minL')

                # next loop over contractions
                for basisContraction in basisGroup.xpath('molpro-output:basisContraction',namespaces=namespaces):
                    kntcontr += 1
                    cc=np.float64(re.sub(' +',' ',basisContraction.text.replace('\n','').lstrip().rstrip()).split(" "))
                    nprim = 0
                    for c in cc:
                        if c != 0.:
                            nprim += 1
                    #sys.stdout.write(('cc = '+[str(cci) for cci in cc] + "\n"))   # Prints the contraction coefficients
                    baseatom.append(nprim)
                    baseatom.append(lcontr)
                    for i in range((2*int(lcontr)+1)):
                        lvec.append(int(lcontr))
                    for i in range(len(alpha)):
                        if cc[i] != 0.:
                            baseatom.append(alpha[i])
                    for i in range(len(alpha)):
                        if cc[i] != 0.:
                            baseatom.append(cc[i])
        datosbase.append(kntcontr)
        #sys.stdout.write( "baseatom = " + baseatom + "\n")
        for aux in baseatom:
            datosbase.append(aux)
        iatom+=1
    #sys.stdout.write( "lvec = " + lvec + "\n")
    #sys.stdout.write( "datosbase = " + datosbase + "\n")
    knt = 0
    facts = [0.5*sqrt(math.acos(-1))]
    for i in range(1,20):
        facts.append(facts[-1] * 0.5 * (2*i+1))
    #sys.stdout.write( "facts = " + facts + "\n")
    while knt < len(datosbase):
        numcontr = datosbase[knt]
        #sys.stdout.write( "\n" + numcontr + "\n")
        escritor.writerow([])
        escritor.writerow([numcontr])
        knt +=1
        for j in range(numcontr):
            nprim = datosbase[knt]
            nbas += (2*int(datosbase[knt+1])+1)
            l = int(datosbase[knt+1])
            #sys.stdout.write( nprim + l + "\n")
            escritor.writerow([nprim, datosbase[knt+1]])
            knt += 2
            aux = [ e for e in datosbase[knt:knt+nprim] ]
            for k in range(len(aux)//5):
                escritor.writerow([("%.15e" %bux) for bux in aux[k*5:(k+1)*5]])
            if (len(aux)%5 != 0):
                escritor.writerow([("%.15e" %bux) for bux in aux[(len(aux)//5)*5:len(aux)+1]])
            renorm = [sqrt(sqrt(4.*math.pow(2.*float(e),(2*l+3)))/facts[l]) for e in datosbase[knt:knt+nprim]]
            #sys.stdout.write( "renorm = " + renorm + "\n")
            #sys.stdout.write((" ".join( str(repr(e)) for e in datosbase[knt:knt+nprim] )) + "\n"))
            knt += nprim
            aux = [ float(datosbase[m])*renorm[m-knt] for m in range(knt,knt+nprim) ]
            for k in range(len(aux)//5):
                escritor.writerow([("%.15e" %bux) for bux in aux[k*5:(k+1)*5]])
            if (len(aux)%5 != 0):
                escritor.writerow([("%.15e" %bux) for bux in aux[(len(aux)//5)*5:len(aux)+1]])
            #sys.stdout.write((" ".join( str(repr(e)) for e in datosbase[knt:knt+nprim] )) + "\n"))
            #escritor.writerow(aux)
            knt += nprim
    nbasis = nbas
    return nbas # end basis evaluation

def writeOrbitals(nbas, molecule): # Prints molecular orbitals
    orbitalSets = molecule.xpath('molpro-output:orbitals',namespaces=namespaces)
    if len(orbitalSets) != 1:
        raise Exception('something is wrong: there should be just one orbital set')
    orbitalSet = orbitalSets[0]
    mos = []
    occv = []
    orderorb = []
    symorb = []
    kntorb = 0
    for orbital in orbitalSet.xpath('molpro-output:orbital',namespaces=namespaces):
        occ = np.float64(orbital.get('occupation'))
        orderorb.append([kntorb,np.float64(orbital.get('energy'))])
        kntorb += 1
        symorb.append(np.float64(orbital.get('state_symmetry')))
        occv.append(occ)
        mosi = array(re.sub(' +',' ',orbital.text.lstrip().rstrip().replace('\n','')).split(" "))
        #sys.stdout.write( 'orb no. ' + kntorb + "\n"))
        #sys.stdout.write( mosi  + "\n"))
        mosround = []
        for bux in mosi:
            if (abs(float(bux)) < 1.e-12):
                bux = 0.
            mosround.append(bux)
        knt = 0
        #sys.stdout.write( "len(lvec) = " + len(lvec) + " len(mosround) = " + len(mosround) + "\n"))
        mosi = []
        while(knt < len(mosround)):
            #sys.stdout.write( "lvec[",knt,"] = " + lvec[knt] + "\n"))
            if (lvec[knt] == 0):
                mosi.append(float(mosround[knt]))
                knt += 1
            elif (lvec[knt] == 1):
                mosi.append(float(mosround[knt+1]))  # loads p_y in p_{-1}
                mosi.append(float(mosround[knt+2]))  # loads p_z in p_{0}
                mosi.append(float(mosround[knt]))    # loads p_x in p_{+1}
                knt += 3
            elif (lvec[knt] == 2):
                mosi.append(float(mosround[knt+1]))  # loads d_{-2}
                mosi.append(float(mosround[knt+4]))  # loads d_{-1}
                mosi.append(float(mosround[knt]))    # loads d_{0}
                mosi.append(float(mosround[knt+2]))  # loads d_{+1}
                mosi.append(float(mosround[knt+3]))  # loads d_{+2}
                knt += 5
            elif (lvec[knt] == 3):
                mosi.append(float(mosround[knt+5]))  # loads f_{-3}
                mosi.append(float(mosround[knt+4]))  # loads f_{-2}
                mosi.append(float(mosround[knt+1]))  # loads f_{-1}
                mosi.append(float(mosround[knt+2]))  # loads f_{0}
                mosi.append(float(mosround[knt]))    # loads f_{+1}
                mosi.append(float(mosround[knt+6]))  # loads f_{+2}
                mosi.append(float(mosround[knt+3]))  # loads f_{+3}
                knt += 7
            elif (lvec[knt] == 4):
                mosi.append(float(mosround[knt+6]))  # loads g_{-4}
                mosi.append(float(mosround[knt+8]))  # loads g_{-3}
                mosi.append(float(mosround[knt+1]))  # loads g_{-2}
                mosi.append(float(mosround[knt+4]))  # loads g_{-1}
                mosi.append(float(mosround[knt]))    # loads g_{0}
                mosi.append(float(mosround[knt+2]))  # loads g_{+1}
                mosi.append(float(mosround[knt+5]))  # loads g_{+2}
                mosi.append(float(mosround[knt+7]))  # loads g_{+3}
                mosi.append(float(mosround[knt+3]))  # loads g_{+4}
                knt += 9
            elif (lvec[knt] == 5):
                mosi.append(float(mosround[knt+7]))  # loads h_{-5}
                mosi.append(float(mosround[knt+4]))  # loads h_{-4}
                mosi.append(float(mosround[knt+5]))  # loads h_{-3}
                mosi.append(float(mosround[knt+10])) # loads h_{-2}
                mosi.append(float(mosround[knt+1]))  # loads h_{-1}
                mosi.append(float(mosround[knt+8]))  # loads h_{0}
                mosi.append(float(mosround[knt]))    # loads h_{+1}
                mosi.append(float(mosround[knt+2]))  # loads h_{+2}
                mosi.append(float(mosround[knt+3]))  # loads h_{+3}
                mosi.append(float(mosround[knt+6]))  # loads h_{+4}
                mosi.append(float(mosround[knt+9]))  # loads h_{+5}
                knt += 11
            elif (lvec[knt] == 6):
                mosi.append(float(mosround[knt+6]))  # loads i_{-6}
                mosi.append(float(mosround[knt+4]))  # loads i_{-5}
                mosi.append(float(mosround[knt+8]))  # loads i_{-4}
                mosi.append(float(mosround[knt+10])) # loads i_{-3}
                mosi.append(float(mosround[knt+1]))  # loads i_{-2}
                mosi.append(float(mosround[knt+11])) # loads i_{-1}
                mosi.append(float(mosround[knt+9]))  # loads i_{0}
                mosi.append(float(mosround[knt+12])) # loads i_{+1}
                mosi.append(float(mosround[knt+5]))  # loads i_{+2}
                mosi.append(float(mosround[knt+7]))  # loads i_{+3}
                mosi.append(float(mosround[knt+3]))  # loads i_{+4}
                mosi.append(float(mosround[knt+2]))  # loads i_{+5}
                mosi.append(float(mosround[knt]))    # loads i_{+6}
                knt += 13
            else:
                sys.stdout.write( "Not prepared for l > 6" + "\n")
                return
        mos.append(mosi)
    #sys.stdout.write( "unordered orderorb " + orderorb     + "\n"))
    orderorb = sorted(orderorb, key=lambda orb: orb[1])
    #sys.stdout.write( "ordered orderorb " + orderorb + "\n"))
    iord = [a[0] for a in orderorb]    # Orbital indices ordered on increasing energy
    #sys.stdout.write( "iord = ", iord
    #sys.stdout.write( "nbas = " , nbas, " knt = ", knt, " len(mosi) = ", len(mosi)
    escritor.writerow([nbas,kntorb,kntorb])
    for i in range(kntorb):
        escritor.writerow([])
        for j in range(nbas//5):
            aux = [mos[iord[i]][k] for k in range(j*5,j*5+5)]
            escritor.writerow([("%.15e" %bux) for bux in aux])
        if (nbas%5 != 0):
#            sys.stdout.write( "(nbas//5)*5 = " + str((nbas//5)*5) + "  nbas = " + str(nbas) + "\n")
#            sys.stdout.write( "k = ")
#            [sys.stdout.write(str(k) + " ") for k in range((nbas//5)*5,nbas)]
            aux = [mos[iord[i]][k] for k in range((nbas//5)*5,nbas)]
#            sys.stdout.write( "aux = " + str(aux) + "\n")
            escritor.writerow([("%.15e" %bux) for bux in aux])
        escritor.writerow([])
    #sys.stdout.write( "kntorb = " + kntorb + "\n"))
    den = []
    escritor2.writerow([nbas])
    for i in range(nbas):
        for j in range(i+1):
            aux = 0.
            for k in range(kntorb):
                aux += occv[k]*mos[k][i]*mos[k][j]
                if (abs(aux) < 1.e-16):
                    aux = 0.
            den.append(aux)
    #sys.stdout.write( nbas + len(den) + "\n"))
    for i in range((len(den)//5)):
        escritor2.writerow([("%.15e" %aux) for aux in den[i*5:i*5+5]])
    if (len(den)/5 != 0):
        escritor2.writerow([("%.15e" %aux) for aux in den[(len(den)//5)*5:len(den)]])
    return
        
#  THE MAIN PROGRAM
#  ================


molecules = document.xpath('/*/molpro-output:molecule',namespaces=namespaces)
if len(molecules) == 0:
    fout.write(inpath+importfile+" is not a MOLPRO xml file suitable for DAMQT\n")
    fout.write("Check that it has been generated with the option {put,xml,"+importfile+";keepspherical}\n")
    fout.write("included in the .com input file\n")
    fout.close()
for answer in document.xpath('//molpro-output:variables/molpro-output:variable[@name="_ANGSTROM"]/molpro-output:value',namespaces=namespaces):
    Angstrom=np.float64(answer.text)
    #sys.stdout.write(('Angstrom = ' + Angstrom + "\n")))
    

fout.write("MOLPRO_interface output\n")
fout.write("==========================\n\n")
fout.write("inpath = "+inpath+"\n")
fout.write("outpath = "+outpath+"\n")
fout.write("PROJECT NAME: "+projectname+"\n")
#sys.stdout.write( "molecules = " + molecules + "\n"))
for molecule in molecules:
    for answer in molecule.xpath('//molpro-output:variables/molpro-output:variable[@name="_ANGSTROM"]/molpro-output:value',namespaces=namespaces):
        Angstrom=np.float64(answer.text)
    if len(molecules) > 1:
        file = outpath+projectname+'_'+ molecule.get('id')
        fout.write("\nData for molecule: "+molecule.get('id')+"\n")
        fout.write("--------------------------------\n")
    else:
        file = outpath+projectname
    f=open(file+'.ggbs','w')
    escritor = csv.writer(f, delimiter=' ', quoting=csv.QUOTE_NONE)
    atoms = molecule.xpath('cml:atomArray/cml:atom',namespaces=namespaces)# Kept to provide compatibility with old molpro's xml files
    if len(atoms) == 0:
        atoms = molecule.xpath('cml:molecule/cml:atomArray/cml:atom',namespaces=namespaces)
    if len(atoms) == 0:
        fout.write("Error reading atoms from cml:molecule/cml:atomArray/cml:atom\n")
        raise Exception('Error reading atoms from cml:molecule/cml:atomArray/cml:atom')
    #sys.stdout.write( "atoms = " + atoms + "\n"))
    escritor.writerow([len(atoms)])
    elementTypes=empty((len(atoms)),dtype='a4')
    iatom = 0
    aux = []
    centrocargas = [0,0,0]
    cargatotal = 0
    for atom in atoms:
        aux.append([np.float64(atom.get('x3'))*Angstrom,np.float64(atom.get('y3'))*Angstrom,np.float64(atom.get('z3'))*Angstrom, atomicnumber[atom.get('elementType')]])
        centrocargas = [centrocargas[j]+aux[iatom][j]*atomicnumber[atom.get('elementType')] for j in range(3)]
        cargatotal += atomicnumber[atom.get('elementType')]
        #sys.stdout.write( str(np.float64(atom.get('x3'))*Angstrom)+str(np.float64(atom.get('y3'))*Angstrom)+str(np.float64(atom.get('z3'))*Angstrom) + str(atomicnumber[atom.get('elementType')) + "\n"))]
        iatom+=1
    centrocargas = [0,0,0]   # Keeps original coordinates without translation of origin to center of positive charges
    for i in range(iatom):
        for j  in range(3):
            aux[i][j] = aux[i][j]-centrocargas[j]/cargatotal
    for i in range(iatom):
        escritor.writerow(aux[i])
    numcontr = 0
    nbas = writeBasis(molecule)
    f.close()
    f=open(file+'.GAorba','w')
    escritor = csv.writer(f, delimiter=' ', quoting=csv.QUOTE_NONE)
    f2=open(file+'.den','w')
    escritor2 = csv.writer(f2, delimiter=' ', quoting=csv.QUOTE_NONE)
    writeOrbitals(nbas,molecule)
    f.close()
    f2.close()
    fout.write("Number of centers: "+str(iatom)+"\n")
    fout.write("Number of contractions: "+str(nbas)+"\n")
fout.close()
