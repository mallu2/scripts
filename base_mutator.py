#!/usr/bin/python
# -*- coding: utf-8 -*- 

import __main__
__main__.pymol_argv = [ 'pymol', '-qc']

import os, sys, math, random, argparse, time, shutil
import pymol     
from pymol import cmd

pymol.finish_launching()

__author__ = 'Dennis M. Krueger'
__copyright__ = 'Copyright 2017, Dennis M. Krueger'
__version__ = '1.0'

class MyParser(argparse.ArgumentParser):
  def error(self, message): 
    sys.stderr.write('error: %s\n' % message)
    self.print_help()
    sys.exit(2)

def check_file(x):
  if os.path.isfile(x) == False:               
    raise argparse.ArgumentTypeError("File does not exist") 
  return x

parser=MyParser(description='Syntax description')
parser.add_argument('-p', type=check_file, help='PDB file ', required=True)
parser.add_argument('-s', type=check_file, help='Sequence file ', required=True)

#print """
#  ____                   __  __       _        _
# | __ )  __ _ ___  ___  |  \/  |_   _| |_ __ _| |_ ___  _ __
# |  _ \ / _` / __|/ _ \ | |\/| | | | | __/ _` | __/ _ \| '__|
# | |_) | (_| \__ \  __/ | |  | | |_| | || (_| | || (_) | |
# |____/ \__,_|___/\___| |_|  |_|\__,_|\__\__,_|\__\___/|_|
#
# -. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .
# ||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|
# |/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\||
# ~   `-~ `-`   `-~ `-`   `-~ `-~   `-~ `-`   `-~ `-`   `-~ `-
# ____________________________________________________________\n
#        written by Dennis M. Krueger, September 2017
#         ICM, Uppsala University, Uppsala, Sweden\n   
#"""

try:
  args = parser.parse_args()
except:
  print("\n")
  sys.exit()
  cmd.quit()

sFileS = args.s
sFileP = args.p
sPath = os.getcwd()

sBaseData_DA = [
'ATOM      1  C1\'  DA     1      -1.898 118.540  -7.924  1.00  0.00           C',
'ATOM      2  N1   DA     1      -3.988 117.321  -3.422  1.00  0.00           N',
'ATOM      3  C2   DA     1      -2.735 117.271  -3.888  1.00  0.00           C',
'ATOM      4  N3   DA     1      -2.270 117.599  -5.089  1.00  0.00           N',
'ATOM      5  C4   DA     1      -3.266 118.029  -5.883  1.00  0.00           C',
'ATOM      6  C5   DA     1      -4.605 118.137  -5.553  1.00  0.00           C',
'ATOM      7  C6   DA     1      -4.969 117.759  -4.241  1.00  0.00           C',
'ATOM      8  N6   DA     1      -6.216 117.809  -3.762  1.00  0.00           N',
'ATOM      9  N7   DA     1      -5.336 118.622  -6.631  1.00  0.00           N',
'ATOM     10  C8   DA     1      -4.435 118.795  -7.572  1.00  0.00           C',
'ATOM     11  N9   DA     1      -3.159 118.463  -7.183  1.00  0.00           N']

sBaseData_DG = [
'ATOM      1  C1\'  DG     1      14.428  -7.039  -9.545  1.00  0.00           C',
'ATOM      2  N1   DG     1      10.104  -5.997 -12.002  1.00  0.00           N',
'ATOM      3  C2   DG     1      11.061  -5.135 -11.514  1.00  0.00           C',
'ATOM      4  N2   DG     1      10.857  -3.825 -11.736  1.00  0.00           N',
'ATOM      5  N3   DG     1      12.142  -5.527 -10.859  1.00  0.00           N',
'ATOM      6  C4   DG     1      12.192  -6.868 -10.721  1.00  0.00           C',
'ATOM      7  C5   DG     1      11.290  -7.816 -11.166  1.00  0.00           C',
'ATOM      8  C6   DG     1      10.135  -7.382 -11.877  1.00  0.00           C',
'ATOM      9  O6   DG     1       9.222  -8.070 -12.367  1.00  0.00           O',
'ATOM     10  N7   DG     1      11.704  -9.093 -10.816  1.00  0.00           N',
'ATOM     11  C8   DG     1      12.832  -8.905 -10.182  1.00  0.00           C',
'ATOM     12  N9   DG     1      13.196  -7.581 -10.110  1.00  0.00           N']

sBaseData_DT = [
'ATOM      1  C1\'  DT     1      -5.939 126.836  -5.838  1.00  0.00           C',
'ATOM      2  N1   DT     1      -5.790 125.997  -4.632  1.00  0.00           N',
'ATOM      3  C2   DT     1      -4.624 125.275  -4.511  1.00  0.00           C',
'ATOM      4  O2   DT     1      -3.752 125.272  -5.368  1.00  0.00           O',
'ATOM      5  N3   DT     1      -4.513 124.558  -3.345  1.00  0.00           N',
'ATOM      6  C4   DT     1      -5.428 124.490  -2.311  1.00  0.00           C',
'ATOM      7  O4   DT     1      -5.174 123.815  -1.312  1.00  0.00           O',
'ATOM      8  C5   DT     1      -6.640 125.255  -2.513  1.00  0.00           C',
'ATOM      9  C6   DT     1      -6.761 125.956  -3.651  1.00  0.00           C',
'ATOM     10  C7   DT     1      -7.695 125.235  -1.453  1.00  0.00           C']

sBaseData_DC = [
'ATOM      1  C1\'  DC     1      -1.323 109.699  -3.137  1.0  0.00           C',
'ATOM      2  N1   DC     1      -2.591 109.952  -3.837  1.00  0.00           N',
'ATOM      3  C2   DC     1      -3.612 110.589  -3.137  1.00  0.00           C',
'ATOM      4  O2   DC     1      -3.423 110.887  -1.947  1.00  0.00           O',
'ATOM      5  N3   DC     1      -4.775 110.869  -3.769  1.00  0.00           N',
'ATOM      6  C4   DC     1      -4.941 110.527  -5.046  1.00  0.00           C',
'ATOM      7  N4   DC     1      -6.106 110.840  -5.624  1.00  0.00           N',
'ATOM      8  C5   DC     1      -3.922 109.855  -5.784  1.00  0.00           C',
'ATOM      9  C6   DC     1      -2.772 109.589  -5.145  1.00  0.00           C']

fDA = open("DA.pdb","w")
for i in range(len(sBaseData_DA)):
  fDA.write(sBaseData_DA[i]+"\n")
fDA.close()

fDG = open("DG.pdb","w")
for i in range(len(sBaseData_DG)):
  fDG.write(sBaseData_DG[i]+"\n")
fDG.close()

fDT = open("DT.pdb","w")
for i in range(len(sBaseData_DT)):
  fDT.write(sBaseData_DT[i]+"\n")
fDT.close()

fDC = open("DC.pdb","w")
for i in range(len(sBaseData_DC)):
  fDC.write(sBaseData_DC[i]+"\n")
fDC.close()

aAtoms = ["P","OP1","OP2","O3'","O5'","C1'","C2'","C3'","C4'","C5'","O3'","O4'","O5'","N1","C6","N9","C8"]
newPDB = []

sSeqS = ""
sSeqP = ""
iSeqL = 0

f1 = open(sFileS,"r")
sSeqS += f1.readline().strip()
iSeqL = len(sSeqS)
sSeqS += f1.readline().strip()[::-1]
f1.close()
sSeqS = sSeqS.replace(" ","").upper()

asSeqS = list(sSeqS)
for i in range(len(asSeqS)):
  if(asSeqS[i] not in ["A","T","G","C"]):
    print(" ERROR: Wrong base identifier. "+asSeqS[i])
    sys.exit()
    cmd.quit()

aNoBases = []

aID = []
f2 = open(sFileP,"r")
for line in f2.readlines():
  baseID = line[23:26].strip()
  atomName = line[13:17].strip()
  baseName = line[17:20].strip()[:2]
  chain = line[21:22]
  
  if(line[:4] == "ATOM"):
    if(baseName in ["DA","DT","DG","DC"]):
      aID.append(baseID+" "+chain)
  
      if(atomName in aAtoms):
        if(atomName == "N9" or atomName == "C8"):
          if(baseName == "DA" or baseName == "DG"):
            newPDB.append(line)
        if(atomName == "N1" or atomName == "C6"):  
          if(baseName == "DC" or baseName == "DT"):
            newPDB.append(line)
        if(atomName != "N9" and atomName != "C8" and atomName != "N1" and atomName != "C6"):
          newPDB.append(line)
  if((line[:4] == "ATOM" or line[:6] == "HETATM" or line[:3] == "TER") and baseName not in ["DA","DT","DG","DC"]):
    aNoBases.append(line)
f2.close()

aID = list(set(aID))

if(len(aID) != len(sSeqS)):
  print(" ERROR: Sequences have different lengths. PDB file:"+str(len(aID))+" Sequence file:"+str(len(sSeqS)))
  sys.exit()
  cmd.quit()

fnewPDB = open(sFileP[:-4]+"_mutated.pdb","w")
for i in range(len(newPDB)):
  fnewPDB.write(newPDB[i])
fnewPDB.close()

bB = True
baseID = ""
iZ = -1

fFinal = open(sFileP[:-4]+"_mod.pdb","w")

for i in range(len(newPDB)):
  
  if(bB == True):  
    fTmp = open("tmp_vec.pdb","w")
    baseID = newPDB[i][23:26].strip()  
    chainID = line[21:22]
    bB = False
  
  if(newPDB[i][23:26].strip() == baseID):
    fTmp.write(newPDB[i])
    sBasePDB = newPDB[i][17:20].strip()
    sChainPDB = newPDB[i][21:22].strip()
  
  if((newPDB[i][23:26].strip() != baseID and newPDB[i][17:20].strip() != chainID) or i == (len(newPDB)-1)): 
    fTmp.close()

    iZ +=1
    sBaseSeq = "D"+asSeqS[iZ]
    sFileName = sBaseSeq+".pdb"
    
    cmd.feedback("disable","all","actions")
    cmd.feedback("disable","all","results")
    cmd.load(sFileName,"base")
    cmd.load("tmp_vec.pdb","mutated")
    
    if(sBaseSeq == "DA" or sBaseSeq == "DG"):
      sele1 = ["C1'","N9", "C8"]
    else:
      sele1 = ["C1'","N1", "C6"]
                   
    if(sBasePDB == "DA" or sBasePDB == "DG" or sBasePDB == "DG5"):
      sele2 = ["C1'","N9", "C8"]
    else:
      sele2 = ["C1'","N1", "C6"]
    
    cmd.pair_fit("base///"+sBaseSeq+"/"+sele1[0],"mutated///"+sBasePDB+"/"+sele2[0],"base///"+sBaseSeq+"/"+sele1[1],"mutated///"+sBasePDB+"/"+sele2[1],"base///"+sBaseSeq+"/"+sele1[2],"mutated///"+sBasePDB+"/"+sele2[2])    
        
    cmd.save("tmp_base.pdb","base")
    cmd.reinitialize()
    
    aOP = []
        
    ff = open("tmp_vec.pdb","r")
    for line in ff.readlines():
      if(line[:4] == "ATOM"):
        atomName = line[13:17].strip()
        if(atomName not in ["C1'","N9","N1","C6","C8"]):       

          sZ = str(iZ+1)
          if(iZ+1 < 10):
            sZ = "    "+sZ
          if(iZ+1 >= 10 and iZ+1 < 100):
            sZ = "   "+sZ
          if(len(sBaseSeq) == 2):
            sBaseSeq = " "+sBaseSeq
          
          if(atomName not in ["OP1","OP2"]):
            fFinal.write(line[:17]+sBaseSeq+" "+sZ+line[27:])
          else:
            aOP.append(line[:17]+sBaseSeq+" "+sZ+line[27:])
    ff.close()
    
    ff = open("tmp_base.pdb","r")     
    for line in ff.readlines():
      if(line[:4] == "ATOM"):
          atomName = line[13:17].strip()
          sZ = str(iZ+1)
          if(iZ+1 < 10):
            sZ = "    "+sZ
          if(iZ+1 >= 10 and iZ+1 < 100):
            sZ = "   "+sZ
          fFinal.write(line[:21]+sZ+line[27:])
    
    if(len(aOP) > 0):
      fFinal.write(aOP[0])
      fFinal.write(aOP[1])
    ff.close()
        
    if(i != (len(newPDB)-1)):
      fTmp = open("tmp_vec.pdb","w") 
      baseID = newPDB[i][23:26].strip()
      fTmp.write(newPDB[i])
      
    if(iZ+1 == iSeqL):
      fFinal.write("TER\n") 

fFinal.write("TER\n")

for i in range(len(aNoBases)):
  fFinal.write(aNoBases[i])
  
fFinal.close()

os.remove(sFileP[:-4]+"_mutated.pdb")
os.remove("tmp_base.pdb")
os.remove("tmp_vec.pdb")
os.remove("DA.pdb")
os.remove("DT.pdb")
os.remove("DG.pdb")
os.remove("DC.pdb")

print("\n Bases mutated successfully!\n")
