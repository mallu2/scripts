#/usr/bin/python

import sys
import os
import numpy as np
import re

def RidPlus(file,integer):
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
    for line in file.readlines():
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
    RI = [(float(ri) + integer) for ri in Residueid]    
    X = [float(x) for x in Xcord]    
    Y = [float(y) for y in Ycord]    
    Z = [float(z) for z in Zcord]    

    return(Atom,AI,Atomname,Residuename,RI,X,Y,Z,ele)

def PDB_A(line):
    return('{dl[0]:6s}{dl[1]:5.0f} {dl[2]:4s} {dl[3]:3s} {otf[2]:1s}{dl[4]:4.0f}    {dl[5]:8.3f}{dl[6]:8.3f}{dl[7]:8.3f}{otf[0]:6.2f}{otf[1]:6.2f}          {dl[8]:3s}\n'.format(dl=line,otf=(1,0,'A')))

def PDB_B(line):
    return('{dl[0]:6s}{dl[1]:5.0f} {dl[2]:4s} {dl[3]:3s} {otf[2]:1s}{dl[4]:4.0f}    {dl[5]:8.3f}{dl[6]:8.3f}{dl[7]:8.3f}{otf[0]:6.2f}{otf[1]:6.2f}          {dl[8]:3s}\n'.format(dl=line,otf=(1,0,'B')))

def PDB_C(line):
    return('{dl[0]:6s}{dl[1]:5.0f} {dl[2]:4s} {dl[3]:3s} {otf[2]:1s}{dl[4]:4.0f}    {dl[5]:8.3f}{dl[6]:8.3f}{dl[7]:8.3f}{otf[0]:6.2f}{otf[1]:6.2f}          {dl[8]:3s}\n'.format(dl=line,otf=(1,0,'C')))

def PDB_D(line):
    return('{dl[0]:6s}{dl[1]:5.0f} {dl[2]:4s} {dl[3]:3s} {otf[2]:1s}{dl[4]:4.0f}    {dl[5]:8.3f}{dl[6]:8.3f}{dl[7]:8.3f}{otf[0]:6.2f}{otf[1]:6.2f}          {dl[8]:3s}\n'.format(dl=line,otf=(1,0,'D')))

def find_ends_p(file):
    "Find the atom ID of the respective ends in the chain."
    end = []
    for line in file.readlines():
        if ("OXT") in line:
            end.append(int(line[6:11]))
    return(end)        

def find_ends_d(file):
    "Find the atom ID of the respective ends in the chain."
    end = []
    for line in file.readlines():
        if ("OP2 D") in line:
            end.append(int(line[6:11]))
    return(end)        

"Open PDB file"
f1 = open("{}".format(sys.argv[1]), 'r')
f2 = open("{}".format(sys.argv[2]), 'w')
print('opening file '+f2.name)

"Get the indexes of ends"
ENDS = find_ends_p(f1)
print(ENDS)
f1.close
f1 = open("{}".format(sys.argv[1]), 'r')
ENDSD = find_ends_d(f1)
print(ENDSD)

f1.close()

f1 = open("{}".format(sys.argv[1]), 'r')

"Write out the new data with changed atom numbers"
listA = RidPlus(f1,1)
f1.close()

f1 = open("{}".format(sys.argv[1]), 'r')

listB = RidPlus(f1,(-334))
f1.close()

f1 = open("{}".format(sys.argv[1]), 'r')

"Write out the new data with changed atom numbers for DNA"
listC = RidPlus(f1,(-670))
f1.close()
f1 = open("{}".format(sys.argv[1]), 'r')
listD = RidPlus(f1,(-694))
f1.close()

print(len(listA))
print(len(listA[0]))
print(len(listA))
print(len(listA[0]))

"Write new chain A"
for i in list(range(0,ENDS[0])):
    line = [listA[0][i],listA[1][i],listA[2][i],listA[3][i],listA[4][i],listA[5][i],listA[6][i],listA[7][i],listA[8][i]]
    f2.write(PDB_A(line))
"Write new chain B"
for i in range(ENDS[0],ENDS[1]):
    line = [listB[0][i],listB[1][i],listB[2][i],listB[3][i],listB[4][i],listB[5][i],listB[6][i],listB[7][i],listB[8][i]]
    f2.write(PDB_B(line))
"Write new chain C"
for i in range(ENDS[1],ENDSD[0]):
    line = [listC[0][i],listC[1][i],listC[2][i],listC[3][i],listC[4][i],listC[5][i],listC[6][i],listC[7][i],listC[8][i]]
    f2.write(PDB_C(line))
"Write new chain D"
for i in range(ENDSD[0],ENDSD[1]):
    line = [listD[0][i],listD[1][i],listD[2][i],listD[3][i],listD[4][i],listD[5][i],listD[6][i],listD[7][i],listD[8][i]]
    f2.write(PDB_D(line))

f1.close()
f2.close()

"Replace each terminus with a TER"
 
def TER_ends(file1,file2):
    "Find the atom ID of the respective ends in the chain."
    end = []
    for line in file1.readlines():
        if 'OXT' in line:
            line = 'TER\n'
        file2.write(line)

f2 = open("{}".format(sys.argv[2]), 'r')
f3 = open("{}".format(sys.argv[3]), 'w')
TER_ends(f2,f3)

f2.close()
f3.close()



