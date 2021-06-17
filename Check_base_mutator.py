#/usr/bin/python

import sys
import os
import numpy as np
import re

def RidPlus(file,start,end):
    """Rertuns the residue-id modified by integer."""
    Atom = []
    Atomid = []
    Atomname = []
    Residuename = []
    Residueid = []
    Xcord = []
    Ycord = []
    Zcord = []
    ele = []
    print(file)
    lines=file.readlines()
    #print(lines)
    for line in lines[start:end]:
        print(line)
        Atom.append(line[0:6])
        Atomid.append(line[6:11])
        Atomname.append(line[12:16])
        Residuename.append(line[17:20])
        Residueid.append(line[22:26])
        Xcord.append(line[30:38])
        Ycord.append(line[38:46])
        Zcord.append(line[46:54])
        ele.append(line[76:80])
    AI = [float(a) for a in Atomid]    
    RI = [float(ri) for ri in Residueid]    
    X = [float(x) for x in Xcord]    
    Y = [float(y) for y in Ycord]    
    Z = [float(z) for z in Zcord]    
    #print(Z)
    return(Atom,AI,Atomname,Residuename,RI,X,Y,Z,ele)

def PDB_A(line):
    return('{dl[0]:6s}{dl[1]:5.0f} {dl[2]:4s} {dl[3]:3s} {otf[2]:1s}{dl[4]:4.0f}    {dl[5]:8.3f}{dl[6]:8.3f}{dl[7]:8.3f}{otf[0]:6.2f}{otf[1]:6.2f}          {dl[8]:3s}'.format(dl=line,otf=(1,0,'A')))

def PDB_B(line):
    return('{dl[0]:6s}{dl[1]:5.0f} {dl[2]:4s} {dl[3]:3s} {otf[2]:1s}{dl[4]:4.0f}    {dl[5]:8.3f}{dl[6]:8.3f}{dl[7]:8.3f}{otf[0]:6.2f}{otf[1]:6.2f}          {dl[8]:3s}'.format(dl=line,otf=(1,0,'B')))

def PDB_C(line):
    return('{dl[0]:6s}{dl[1]:5.0f} {dl[2]:4s} {dl[3]:3s} {otf[2]:1s}{dl[4]:4.0f}    {dl[5]:8.3f}{dl[6]:8.3f}{dl[7]:8.3f}{otf[0]:6.2f}{otf[1]:6.2f}          {dl[8]:3s}'.format(dl=line,otf=(1,0,'C')))

def PDB_D(line):
    return('{dl[0]:6s}{dl[1]:5.0f} {dl[2]:4s} {dl[3]:3s} {otf[2]:1s}{dl[4]:4.0f}    {dl[5]:8.3f}{dl[6]:8.3f}{dl[7]:8.3f}{otf[0]:6.2f}{otf[1]:6.2f}          {dl[8]:3s}'.format(dl=line,otf=(1,0,'D')))

def find_ENDS(file):
    "Find the atom ID of the respective ends in the chain."
    END = []
    lines=file.readlines()
    #print(file)
    Nlines=len(lines)
    print(Nlines)
    for i in range(Nlines):
        line=lines[i]
        if ("MODEL  ") in line:
            END.append(i+1)
        if ("929  O   SER") in line:
            END.append(i+1)
        if ("1858  O   SER") in line:
            END.append(i+1)
        if ("2613  H3T DC3") in line:
            END.append(i+1)
        if ("3380  H3T DG3") in line:
            END.append(i+1)
    ##print(END)
    return(END)        

"Open PDB file"
f1 = open("{}".format(sys.argv[1]), 'r')
f2 = open("{}".format(sys.argv[2]), 'w')
print('opening file '+f2.name)
print('opening file '+f1.name)

ENDS=find_ENDS(f1)
f1.close()
models=sys.argv[3]
f2.write('CRYST1  130.924  130.924  130.924  60.00  60.00  90.00 P 1 \n')
for ii in range(int(models)):
    f1 = open("{}".format(sys.argv[1]), 'r')
    "Write new chain A"
    f2.write('MODEL      {:3}\n'.format(ii+1))
    listA=RidPlus(f1,ENDS[ii*4+ii],ENDS[ii*4+ii+1])
    f1.close()
    for i in range(len(listA[0])):
        line = [listA[0][i],listA[1][i],listA[2][i],listA[3][i],listA[4][i],listA[5][i],listA[6][i],listA[7][i],listA[8][i]]
        f2.write(PDB_A(line))
    f2.write('TER\n')
    
    "Write new chain B"
    f1 = open("{}".format(sys.argv[1]), 'r')
    listB=RidPlus(f1,ENDS[ii*4+ii+1],ENDS[ii*4+ii+2])
    f1.close()
    for i in range(len(listB[0])):
        line = [listB[0][i],listB[1][i],listB[2][i],listB[3][i],listB[4][i],listB[5][i],listB[6][i],listB[7][i],listB[8][i]]
        f2.write(PDB_B(line))
    f2.write('TER\n')
    
    "Write new chain C"
    f1 = open("{}".format(sys.argv[1]), 'r')
    listC=RidPlus(f1,ENDS[ii*4+ii+2],ENDS[ii*4+ii+3])
    f1.close()
    for i in range(len(listC[0])):
        line = [listC[0][i],listC[1][i],listC[2][i],listC[3][i],listC[4][i],listC[5][i],listC[6][i],listC[7][i],listC[8][i]]
        f2.write(PDB_C(line))
    f2.write('TER\n')
    
    "Write new chain D"
    f1 = open("{}".format(sys.argv[1]), 'r')
    listD=RidPlus(f1,ENDS[ii*4+ii+3],ENDS[ii*4+ii+4])
    f1.close()
    for i in range(len(listD[0])):
        line = [listD[0][i],listD[1][i],listD[2][i],listD[3][i],listD[4][i],listD[5][i],listD[6][i],listD[7][i],listD[8][i]]
        f2.write(PDB_D(line))
    f2.write('TER\n') 
    f2.write('ENDMDL\n') 

f1.close()
f2.close()

