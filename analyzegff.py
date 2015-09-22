from __future__ import division
from collections import Counter
from itertools import izip as zip, count
from os import system
from scipy import stats
import numpy as np
import linecache
import re
import argparse
import time
from Bio import SeqIO
import os

# The main function starts
time_start = time.time()

parser = argparse.ArgumentParser(
     prog='analyzegff',
     usage='''python analyzegff.py --gff [Combined gff file] --length [Length of genome] --path [Path of input files] --name [Name of genome]''',
     description='''This program analyzes the gff files from the de novo and repbase runs of repeatmasker''',
     epilog='''It requires numpy, scipy and biopython libraries''')
parser.add_argument('--gff', type=str, help='The combined denovo and repbase gff file', required=True)
parser.add_argument('--length', type=int, help='Length of genome', required=False)
parser.add_argument('--path', type=str, help='Path to input files', required=False)
parser.add_argument('--name', type=str, help='Name of genome', required=False)

args = parser.parse_args()

filename1=args.gff
length2=args.length
path=args.path
name=args.name

numberlines=0

if path == None:
      filename9=filename1
else:
   if path[len(path)-1]=="/":
      filename9=path + filename1
   else:
      filename9=path + "/" + filename1

if path == None:
       filename3=open(filename1,"r")
else:
       filename3=open(os.path.join(path, filename1), "r")      
      
with filename3 as f:
      numberlines=len(f.readlines())
f.close()
length=0
common=0
repbase=0
rep=0
denovo=0
unknown=0
total=0
gypsy=0
copia=0
enspm=0
L1=0
R1=0
hAT=0
LINE=0
SINE=0
DNA=0
simplerepeat=0
helitron=0
retrotransposon=0
mudr=0
penelope=0
satellite=0
rtex=0
unk=0
classified=[]
species=[]
elements=[]

for i in range(numberlines):
   line=linecache.getline(filename9,i+1)    
   line2=line.split()
   if len(line2)!=0:
      start=int(line2[3])
      end=int(line2[4])
      diff=abs(end-start)
      length+=diff
   if "Unknown" in line:
      unknown+=1
      unk+=diff
   if "Gypsy" in line:
      gypsy+=diff
   if "DNA" in line:
      DNA+=diff
   if "Copia" in line:
      copia+=diff
   if "EnSpm" in line:
      enspm+=diff
   if "[" in line:
      rep+=diff
   if "L1" in line:
      L1+=diff
   if "LINE" in line:
      LINE+=diff
   if "SINE" in line:
      SINE+=diff
   if "R1" in line:
      R1+=diff
   if "hAT" in line:
      hAT+=diff
   if "RTE-X" in line:
      rtex+=diff
   if "Simple_repeat" in line:
      simplerepeat+=diff
   if "Helitron" in line:
      helitron+=diff
   if "MuDR" in line:
      mudr+=diff
   if "Penelope" in line:
      penelope+=diff
   if "Satellite" in line:
      satellite+=diff
   if ("LTR Retrotransposon" in line) or ("(LTR)" in line):
     if "Non-LTR Retrotransposon" not in line:
       retrotransposon+=diff
   if "rnd" in line2[8]: # Denovo
      denovo+=1
      abc=line2[8]
      m=line.find('(')
      n=line.find(')')
      elements.append(abc[:m-1])
      tempclass=line[m+1:n]
      if "/" in tempclass:
         q=tempclass.find('/')
         r=len(tempclass)
         classified.append(tempclass[q+1:r])
      else:
         classified.append(tempclass)
   else: # Repbase
      repbase+=1
      abc=line2[8]
      m=line.find('(')
      n=line.find(')')
      elements.append(abc[:m-1])
      o=line.find('[')
      p=line.find(']')
      if o !=-1:
        species.append(line[o+1:p])
      tempclass=line[m+1:n]
      if "/" in tempclass:
         q=tempclass.find('/')
         r=len(tempclass)
         classified.append(tempclass[q+1:r])
      else:
         classified.append(tempclass)
   if ("rnd" in line) and ("[" in line): # Common
      common+=1
      denovo+=1
      repbase+=1
      o=line.find('[')
      p=line.find(']')
      species.append(line[o+1:p])
      abc=line2[8]
      m=line.find('(')
      n=line.find(')')
      elements.append(abc[:m-1])
      tempclass=line[m+1:n]
      if "/" in tempclass:
         q=tempclass.find('/')
         r=len(tempclass)
         classified.append(tempclass[q+1:r])
      else:
         classified.append(tempclass)
      efg=abc[n+2:]
      q=efg.find('(')
      r=efg.find(')')
      elements.append(efg[:q-1])
      tempclass=efg[q+1:r]
      if "/" in tempclass:
         q=tempclass.find('/')
         r=len(tempclass)
         classified.append(tempclass[q+1:r])
      else:
         classified.append(tempclass)
   total+=1
filename3.close()
uniquerepeats=set(classified)
uniquerepeats=list(uniquerepeats)
uniquerepeats=sorted(uniquerepeats)
print uniquerepeats
countrepeats=[[] for i in range(len(uniquerepeats))]
counts1=Counter(classified)
print "Total no. of elements",len(elements)
print "Number of lines",numberlines
uniqueelements=set(elements)
uniqueelements=list(uniqueelements)
uniqueelements=sorted(uniqueelements)
countelements=[[] for i in range(len(uniqueelements))]
counts2=Counter(elements)

i=0
for item in uniquerepeats:
   countrepeats[i]=[str(item),float((int(counts1[item])*100)/len(classified))]
   i+=1 

i=0
for item in uniqueelements:
   countelements[i]=[str(item),float((int(counts2[item])*100)/len(elements)),int(counts2[item])]
   i+=1 

#countrepeats=np.array(countrepeats,dtype=(str,float))
countrepeats=np.array(countrepeats)
countrepeats=countrepeats[np.argsort(countrepeats[:, 1].astype(np.float32))]
countrepeats=countrepeats[::-1] # Reversing array
print countrepeats

countelements=np.array(countelements)
print countelements[:, 1].astype(np.float32)
countelements=countelements[np.argsort(countelements[:, 1].astype(np.float32))]
countelements=countelements[::-1] # Reversing array

uniquespecies=set(species)
uniquespecies=list(uniquespecies)
uniquespecies=sorted(uniquespecies)
countspecies=[[] for i in range(len(uniquespecies))]
counts2=Counter(species)
i=0
for item in uniquespecies:
   countspecies[i]=[str(item),float((int(counts2[item])*100)/len(species))]
   i+=1 

countspecies=np.array(countspecies)
try: 
  countspecies=countspecies[np.argsort(countspecies[:, 1].astype(np.float32))]
  countspecies=countspecies[::-1]   # Reversing array
  print countspecies
except IndexError:
  pass


# Writing the frequency-distribution of different repeats 
if name == None: 
  names="final_output2_" + filename1 + ".txt"
else:
  names="final_output_" + name + ".txt"

with open(names,'w') as xyz:
  xyz.write("Frequency of occurence of different types of repeats:\n")
  for i in range(len(countrepeats)):
    xyz.write("%s\t%s\n" %(countrepeats[i,0],countrepeats[i,1])) 
  
  xyz.write("\n")
# Writing the frequency-distribution of the organismal origin of reference repbase repeats.
  xyz.write("Frequency-distribution of the organismal origin of reference repbase repeats\n")
  for i in range(len(countspecies)):
    xyz.write("%s\t%s\n" %(countspecies[i,0],countspecies[i,1]))  
  xyz.write("\n")
    
  if length2==None:
     xyz.write("%s\n" %("Lengths:"))
     xyz.write("%s\t%i\n" %("Gypsy",gypsy))
     xyz.write("%s\t%i\n" %("Copia",copia))
     xyz.write("%s\t%i\n" %("LTR Retrotransposon",retrotransposon))
     xyz.write("%s\t%i\n" %("Gypsy",gypsy))
     xyz.write("%s\t%i\n" %("Copia",copia))
     xyz.write("%s\t%i\n" %("L1",L1))
     xyz.write("%s\t%i\n" %("R1",R1))
     xyz.write("%s\t%i\n" %("Unknown",unk))
     xyz.write("%s\t%i\n" %("DNA",DNA))
     xyz.write("%s\t%i\n" %("repbase",rep))
     xyz.write("%s\t%i\n" %("hAT",hAT))
     xyz.write("%s\t%i\n" %("LINE",LINE))
     xyz.write("%s\t%i\n" %("SINE",SINE))
     xyz.write("%s\t%i\n" %("Simple repeat",simplerepeat))
     xyz.write("%s\t%i\n" %("Helitron",helitron))
     xyz.write("%s\t%i\n" %("MuDR",mudr))
     xyz.write("%s\t%i\n" %("Penelope",penelope))
     xyz.write("%s\t%i\n" %("Satellite",satellite))
     xyz.write("%s\t%i\n" %("RTE-X",rtex))
  else:
     xyz.write("%s\n" %("Percentages of each type of repeats in genome:"))
     xyz.write("%s\t%f\n" %("Gypsy",(gypsy*100)/length2))
     xyz.write("%s\t%f\n" %("Copia",(copia*100)/length2))
     xyz.write("%s\t%f\n" %("LTR Retrotransposon",(retrotransposon*100)/length2))
     xyz.write("%s\t%f\n" %("L1",(L1*100)/length2))
     xyz.write("%s\t%f\n" %("R1",(R1*100)/length2))
     xyz.write("%s\t%f\n" %("hAT",(hAT*100)/length2))
     xyz.write("%s\t%f\n" %("Unknown",(unk*100)/length2))
     xyz.write("%s\t%f\n" %("repbase",(rep*100)/length2))
     xyz.write("%s\t%f\n" %("LINE",(LINE*100)/length2))
     xyz.write("%s\t%f\n" %("DNA",(DNA*100)/length2))
     xyz.write("%s\t%f\n" %("SINE",(SINE*100)/length2))
     xyz.write("%s\t%f\n" %("Simple repeat",(simplerepeat*100)/length2))
     xyz.write("%s\t%f\n" %("Helitron",(helitron*100)/length2))
     xyz.write("%s\t%f\n" %("MuDR",(mudr*100)/length2))
     xyz.write("%s\t%f\n" %("Penelope",(penelope*100)/length2))
     xyz.write("%s\t%f\n" %("Satellite",(satellite*100)/length2))
     xyz.write("%s\t%f\n" %("RTE-X",(rtex*100)/length2))
  xyz.write("\n")
  xyz.write("Statistics:\n")
  xyz.write("%s %f\n" %("Percentage of unknown denovo hits:",(unknown*100)/total))
  xyz.write("%s %f\n" %("Percentage of common denovo and repbase hits:",(common*100)/total))
  xyz.write("%s %f\n" %("Percentage of denovo hits:",(denovo*100)/total))
  xyz.write("%s %f\n" %("Percentage of repbase hits:",(repbase*100)/total))
  xyz.write("%s %i\n" %("Total length of repeats:",length))
  if length2 != None:
      xyz.write("%s %f\n" %("Percentage of repeats in genome:",(length*100)/length2))        
  xyz.write("%s\t%s\t%s" %("Element","Length","Percentage in genome (printed only if length of genome entered)"))
  xyz.write("\n")
  try:
    blank=uniqueelements.index('')
    del uniqueelements[blank]
  except ValueError:
    pass    
  # Calculating length of elements
  print "Unique elements is:",uniqueelements
  repeatelements=[[] for i in range(len(uniqueelements))]
  s=0
  totallength=0
  #doubles=0
  #check=np.zeros(numberlines)

  #duplicates=np.chararray(numberlines)
  p=0
  for item in uniqueelements: 
     p+=1
     print "Repeat number being analyzed:",p      
     length=0
     elementcount=0
     #linecounter=-1
     relevantelements=[v for v, w in zip(count(), elements) if w == item]
     for i in relevantelements:
         line=linecache.getline(filename9,i+1)
         #linecounter+=1 
         line2=line.split("\t")
         start=int(line2[3])
         end=int(line2[4])
         diff=abs(end-start)
         length+=diff
         elementcount+=1
               #if check[linecounter] > 1:
                    #doubles+=diff
                    #print line,"Item:",item,"No. of times:",check[linecounter],"Previous:",duplicates[linecounter]
               #duplicates[linecounter]=item 
     filename3.close()
     totallength+=length
     repeatelements[s]=[str(item),elementcount,int(length)]
     s+=1
  xyz.write("%s %i\n" %("Total calculated repeat length is:",totallength))
  #xyz.write("%s %i\n" %("Length counted twice is:",doubles))
  #xyz.write("%s %i\n" %("The net repeat length is:",totallength-doubles))
  repeatelements=np.array(repeatelements)
  repeatelements=repeatelements[np.argsort(repeatelements[:, 2].astype(np.int32))]
  repeatelements=repeatelements[::-1]
  for i in range(len(uniqueelements)): # Top elements by genome coverage
     if length2 == None:
         xyz.write("%s\t%i\t%i\n" %(str(repeatelements[i,0]),int(repeatelements[i,1]),int(repeatelements[i,2])))
     else:
         xyz.write("%s\t%i\t%i\t%f\n" %(str(repeatelements[i,0]),int(repeatelements[i,1]),int(repeatelements[i,2]),float((float(repeatelements[i,2])*100)/length2)))

xyz.close()

time_end = time.time()
print "Time taken to run this parser is:",time_end-time_start," seconds"
