# This script requires numpy, matplotlib and biopython libraries
# Script written by Robin Paul

from __future__ import division
from collections import Counter
from os import system
import operator
import csv
from scipy import stats
import numpy as np
import linecache
import re
import argparse
import time
from Bio import SeqIO
import os

time_start = time.time()

# The main function starts

parser = argparse.ArgumentParser(
     prog='comparegff11.py',
     usage='''python comparegff11.py --gff [Combined gff file] --classified [classified fasta file] --repbase [Repbase library file] --path [Path to input files]''',
     description='''This program parses the gff files from the de novo and repbase runs of repeatmasker''',
     epilog='''It requires numpy, matplotlib and biopython libraries''')
parser.add_argument('--gff', type=str, help='The combined denovo and repbase gff file', required=True)
parser.add_argument('--classified', type=str, help='The classified (denovo) fasta file', required=True)
parser.add_argument('--repbase', type=str, help='The repbase library file (e.g "all_plant.ref"', required=True)
parser.add_argument('--path', type=str, help='Path to input files', required=False)


args = parser.parse_args()

filename1=args.gff
filename3=args.classified
filename18=args.repbase
path=args.path
filename4=(filename1[:len(filename1)-8] + ".align")

if path == None:
   filename5=open(filename1,"r")
else:
   filename5=open(os.path.join(path, filename1), "r")

if path == None:
   filename12=open(filename4,"r")
else:
   filename12=open(os.path.join(path, filename4), "r")

with open("all_" + filename1,'w') as abc:
    for row in filename5:
     if row.startswith("##")==True:
        continue
     else:
        abc.write(row)

# getting rid of unwanted "#" out of the alignment file
with open("new_" + filename1 + ".align",'w') as xyz:
  for line in filename12:
        line=line.replace('#','-')
        xyz.write(line)
     
xyz.close()

# Parsing the repbase gff file.

length=0
i=0
detectors=[]

reader=np.loadtxt("all_" + filename1, delimiter='\t',dtype='str')
reader2=np.loadtxt("new_" + filename1 + ".align", delimiter='\t',dtype='str')
interval=np.zeros(len(reader))
alignments=np.zeros(len(reader))
percentidentity=np.zeros(len(reader))

starts=np.zeros(len(reader2)+1)
ends=np.zeros(len(reader2))
starts[0]=0
i=0
k=0
j=0
for line in open("new_" + filename1 + ".align", 'rb'):
    j+=1
    if line.startswith("Gap_init rate")==True:
       starts[i+1]=j
       i+=1
    elif line.startswith("Matrix")==True:
       ends[i]=j
       k+=1
i=0
for row in reader:
       fastanames=row[0]
       start=row[3]
       stop=row[4]
       alignlength=0
       total=0
       strt=starts[i]
       end=ends[i]
       #print "fastanames is",fastanames,start,stop
       for k in range(int(strt)+1,int(end)):
               line=linecache.getline("new_"+ filename1 + ".align", k)
               if (fastanames in line) and (start in line) and (stop in line):
                  #print line
                  continue
               elif len(line) > 27:
                        firstword=line[27:].split(' ')
                        if ("A" in firstword[0]) or ("T" in firstword[0]) or ("G" in firstword[0]) or ("C" in firstword[0]) or ("N" in firstword[0]):
                           alignlength+=len(firstword[0])
                           #print firstword[0]
                        else:
                           iis=line[27:].count("i")
                           dashs=line[27:].count("-")
                           vs=line[27:].count("v")
                           quests=line[27:].count("?")
                           total+=iis+dashs+vs+quests
                           #print iis,dashs,vs,quests
                           #print line[27:]
       alignments[i]=alignlength/2
       try:
          percentidentity[i]=100-(total*200)/alignlength
       except ZeroDivisionError:
          print "row",row,"Alignlength:",alignlength,"total:",total
       s=row[8]
       s2=s[14:]
       s3=s2.find('"')
       detectors.append(s2[:s3])
       diff=int(row[4])-int(row[3]) # Identifying the start and end of the alignment
       interval[i]=diff
       length+=diff
       i+=1


print "MAKE SURE ALL THE NUMBERS IN THE FOLLOWING ARRAY ARE WHOLE NUMBERS!! (CHECK THIS WHEN YOU USE THIS SCRIPT ON A NEW DATASET)"
print alignments

highestidentity=np.zeros(len(alignments))

for i in range(len(alignments)):
    if alignments[i]>=0:
       highestidentity[i]=percentidentity[i]

# Extracting repeats which obey the 80/80 rule

i=-1
unknowns=0
with open("known_" + filename1,'w') as xyz:
 with open("unknown_" + filename1,'w') as pqr:
  for line in open("all_" + filename1,'r'):
    i+=1
    if (alignments[i] >= 0) and (percentidentity[i] >= 0):  # Invoking the 80/80 rule
        xyz.write("%s" %(line)) 
    else:
        pqr.write("%s" %(line))
        unknowns+=1   

with open("denovo_" + filename1,'w') as xyz:
  with open("repbase_" + filename1,'w') as pqr:
    for row in open("known_" + filename1,'r'):
     if row.startswith("##")==True:
        continue
     else:
        row2=row.split()
        s=row2[9]
        s3=s[7:]
        s4=s3.find('"')
        s5=s[7:s4]
        if s5.startswith("rnd") == True:  # Denovo Hit
            xyz.write(row)
        else: # Repbase hit
            row=row.replace('#','-')
            pqr.write(row) 

# Extracting classifications from the (classified) fasta file

if path == None:
   filename15=filename3
else:
   if path[len(path)-1]=="/": 
       filename15=path + filename3
   else:
       filename15=path + "/" + filename3

items=[]
conns=[]
for seq_record in SeqIO.parse(filename15, "fasta"):
     identity=(seq_record.id)
     hashh=identity.find("#")
     imp=identity[hashh+1:len(identity)]
     conn=identity[:hashh]
     conns.append(conn)
     items.append(imp)


reader2=np.loadtxt("repbase_" + filename1, delimiter='\t',dtype='str')
interval=np.zeros(len(reader2))
#print "Length of reader2:", len(reader2)
i=0
detectors2=[]
for row in reader2:
       s=row[8]
       s2=s[14:]
       s3=s2.find('"')
       detectors2.append(s2[:s3])
       i+=1

counts=Counter(detectors2)
unique=set(detectors2)
unique=list(unique)
unique=sorted(unique) # Arranging items in alphabetical order
#print unique
track=np.zeros(len(unique)) # Track which repeats have not been characterized yet           

n=9 # Number of repeats being identified
annot=[]
family=[]
l=-1    

if path == None:
   filename21=filename18
else:
   if path[len(path)-1]=="/": 
       filename21=path + filename18
   else:
       filename21=path + "/" + filename18

if path == None:
   filename22=open(filename18,"r")
else:
   filename22=open(os.path.join(path, filename18), "r")

seq_records1 = SeqIO.parse(filename21, "fasta")
seq_records1 = list(seq_records1)
fastalocations=np.zeros(len(seq_records1))
p=1
q=0
for line in filename22:
    if line.startswith(">") == True:
       fastalocations[q]=p
       q+=1
    p+=1    
        
for item in unique:
    p=0
    for seq_record in seq_records1:
         hits1=re.findall('\\b' + str(item) + '\\b', str(seq_record.id))
         hits2=re.findall('\\b' + str(seq_record.id) + '\\b', str(item))
         if (len(hits1) > 0) and (len(hits2) > 0):
             fastaid=linecache.getline(filename21, int(fastalocations[p]))
             #print item,fastaid
             fastaid2=fastaid
             fastaid=fastaid.split()
             #print fastaid[1],len(fastaid)
             if ("Retrotransposon" in fastaid[0:len(fastaid)]):
                #print fastaid[3]
                annot.append(fastaid[1] + " " + fastaid[2])
                family.append(fastaid[3] + " " + fastaid[4])
                continue
             try:
               annot.append(fastaid[1])
             except IndexError:
               print fastaid
             
             try:
                 family.append(fastaid[2] + " " + fastaid[3] + " " + fastaid[4])
             except IndexError:
                  try:
                       family.append(fastaid[2] + " " + fastaid[3])
                  except IndexError:
                       family.append(fastaid[2])
                       print fastaid2
                 
         p+=1      
             
length=0
i=0
m=0
detectors=[]
print len(unique),len(annot),len(family)
reader1=np.loadtxt("denovo_" + filename1, delimiter='\t',dtype='str')
reader2=np.loadtxt("repbase_" + filename1, delimiter='\t',dtype='str')
interval1=np.zeros(len(reader1))

#denovofastanames=[]
repbasefastanames=[]
#for row1 in reader1:
    #denovofastanames.append(row1[0])

for row2 in reader2:
    repbasefastanames.append(row2[0])
#print repbasefastanames

#counts1=Counter(denovofastanames)
#uniquenames1=set(denovofastanames)

counts2=Counter(repbasefastanames)
uniquenames2=set(repbasefastanames)
uniquenames2=sorted(uniquenames2)
occurrences2=np.zeros(len(uniquenames2))
#print uniquenames2
x=np.zeros(len(uniquenames2))
i=0
for item in uniquenames2:
    x[i]=counts2[item] 
    i+=1
j=1
for k in range(len(x)):
    occurrences2[k]=int(j)
    j+=x[k] 
#for i in range(len(uniquenames2)):
    #print occurrences2[i]
m=0
totals=0
track2=np.zeros(len(reader2)+1,dtype=int) # Finding those repbase results that do not match with denovo results
q=-1
with open("repbasedenovo_" + filename1, 'w') as defg:
 with open("denovo2_" + filename1, 'w') as ijk:
  with open("nestedrepeats_" + filename1, 'w') as pqrs:           
   for row1 in reader1: 
       s1=row1[0]
       start1=int(row1[3])
       stop1=int(row1[4])
       p=0
       abc=0
       for conn in conns:
          if conn in row1[8]:  
               s=row1[8]
               s3=s[14:]
               s4=s3.find('"')
               item2=items[abc]
               connector=conn
          abc+=1
       a=0   
       for item in uniquenames2:
            if (s1 in item) and (item in s1):
              hits=counts2[item]
              nohits=0
              t=-1
              for i in range(int(occurrences2[p]),int(occurrences2[p]+hits)): 
                  row2=linecache.getline("repbase_" + filename1, i)
                  row2=row2.split()
                  t+=1
                  if len(row2)==0:
                       continue
                  start2=int(row2[3])
                  stop2=int(row2[4]) 
                  if ((start2<=start1 and stop1<=stop2) or (start1<=start2 and stop2<=stop1)):
                    l=0
                    pqr=row2[9][7:][:row2[9][7:].find('"')]
                    for item3 in unique:
                        if pqr in item3:
                             repannot=annot[l]
                             repfamily=family[l]
                        l+=1      
                    if abs(start2-stop2)>=abs(start1-stop1):
                        defg.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%s%s%s%s%s%s\n" %(row2[0],row2[1],row2[2],start2,stop2,row2[5],row2[6],row2[7],pqr," (",repannot,") ","[",repfamily,"]"))
                        pqrs.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%s%s%s%s%s%s\n" %(row2[0],row2[1],row2[2],start2,stop2,row2[5],row2[6],row2[7],pqr," (",repannot,") ","[",repfamily,"]"))
                        pqrs.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s %s%s%s\n\n" %(row1[0],row1[1],row1[2],start1,stop1,row1[5],row1[6],row1[7],str(connector),"(",item2,")"))
                    elif abs(start2-stop2)<abs(start1-stop1):
                        defg.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s %s%s%s\n" %(row1[0],row1[1],row1[2],start1,stop1,row1[5],row1[6],row1[7],str(connector),"(",item2,")"))
                        pqrs.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s %s%s%s\n" %(row1[0],row1[1],row1[2],start1,stop1,row1[5],row1[6],row1[7],str(connector),"(",item2,")"))
                        pqrs.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%s%s%s%s%s%s\n\n" %(row2[0],row2[1],row2[2],start2,stop2,row2[5],row2[6],row2[7],pqr," (",repannot,") ","[",repfamily,"]"))  
                        a=1
                    nohits=1
                    m+=1
                    #print "occurrences2[p]:",occurrences2[p],"t:",t
                    #print "occurrences2[p]+t:",occurrences2[p]+t
                    #print "row2:",row2
                    try: 
                          track2[occurrences2[p]+t]=int(1)
                    except IndexError:
                          print occurrences2[p],t
                    totals+=1
                    #print "row1,row2,start1,stop1,start2,stop2",row1[0],row2[0],start1,stop1,start2,stop2
                    
            p+=1
       if a==0:
          ijk.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s %s%s%s\n" %(row1[0],row1[1],row1[2],row1[3],row1[4],row1[5],row1[6],row1[7],str(connector),"(",item2,")"))
defg.close() 
pqrs.close()
with open("repbaseonly_" + filename1,'w') as abcd:
  for i in range(len(reader2)):
      if track2[i]==0:
         row2=linecache.getline("repbase_" + filename1, i)
         if len(row2)==0:
            continue
         else:
            row2=row2.split()
            pqr=row2[9][7:][:row2[9][7:].find('"')]
            l=0
            for item3 in unique:
              if pqr in item3:
                repannot=annot[l]
                repfamily=family[l]
              l+=1 
            abcd.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s %s%s%s%s%s%s\n" %(row2[0],row2[1],row2[2],row2[3],row2[4],row2[5],row2[6],row2[7],pqr,"(",repannot,") ","[",repfamily,"]"))

abcd.close()

filename6="denovo2_" + filename1
filename7="repbasedenovo_" + filename1
filename8="repbaseonly_" + filename1

# Combining denovo2, repbasedenovo and repbaseonly

reader6=np.loadtxt(filename6, delimiter='\t',dtype='str')
reader7=np.loadtxt(filename7, delimiter='\t',dtype='str')
reader8=np.loadtxt(filename8, delimiter='\t',dtype='str')
store6=[]
startnums1=[]
fastanames1=[]
uniquescaffolds1=[]
location1=[]
row4=''
m=-1
for row1 in reader6:
     m+=1
     fastanames1.append(row1[0]) 
     try:
       if str(row1[0])!=str(row4[0]):
          location1.append(m)
          uniquescaffolds1.append(row1[0])
     except IndexError:
          location1.append(m)
          uniquescaffolds1.append(row1[0])
     row4=row1
location1.append(len(reader6))
startnums1=np.zeros(len(uniquescaffolds1)+2)
startnums1[0]=0
i=1
counts1=Counter(fastanames1)
temp=0
for item in uniquescaffolds1:
     startnums1[i]=counts1[item]+temp
     temp=startnums1[i]
     i+=1
startnums1[len(uniquescaffolds1)]=len(fastanames1)   
print "Uniquescaffolds1",uniquescaffolds1[:3]
print "Startnums1",startnums1[:3]

startnums2=[]
store7=[]
fastanames2=[]
uniquescaffolds2=[]
location2=[]
row3=''
m=-1
for row2 in reader7:
     m+=1
     fastanames2.append(row2[0]) 
     #row2=row2.split("\t")
     try:
       if str(row2[0])!=str(row3[0]):
          location2.append(m)
          uniquescaffolds2.append(row2[0])
     except IndexError:
          location2.append(m)
          uniquescaffolds2.append(row2[0])
     row3=row2
location2.append(len(reader7)) 

startnums2=np.zeros(len(uniquescaffolds2)+2)
startnums2[0]=0
i=1
counts2=Counter(fastanames2)
temp=0
for item in uniquescaffolds2:
     startnums2[i]=counts2[item]+temp
     temp=startnums2[i]
     i+=1
startnums2[len(uniquescaffolds2)]=len(fastanames2)  

print "Uniquescaffolds2",uniquescaffolds2[:3]
print "Startnums2",startnums2[:3]

store8=[]
startnums3=[]
fastanames3=[]
uniquescaffolds3=[]
location3=[]
row4=''
m=-1
for row3 in reader8:
     m+=3
     fastanames3.append(row3[0]) 
     try:
       if str(row3[0])!=str(row4[0]):
          location3.append(m)
          uniquescaffolds3.append(row3[0])
     except IndexError:
          location3.append(m)
          uniquescaffolds3.append(row3[0])
     row4=row3
location3.append(len(reader8))

startnums3=np.zeros(len(uniquescaffolds3)+2)
startnums3[0]=0
i=1
counts3=Counter(fastanames3)
temp=0
for item in uniquescaffolds3:
     startnums3[i]=counts3[item]+temp
     temp=startnums3[i]
     i+=1
startnums3[len(uniquescaffolds3)]=len(fastanames3)   
print "Uniquescaffolds3",uniquescaffolds3[:3]
print "Startnums3",startnums3[:3]

fastanames=fastanames1+fastanames2+fastanames3
fastanames=sorted(set(fastanames))

with open(filename1[:len(filename1)-8] + "_final0000.gff",'w') as xyz:
  for item in fastanames:
   
   store6=[]
   for item2 in uniquescaffolds1:
      if (item in item2) and (item2 in item):
         i=uniquescaffolds1.index(item2)
         for j in range(int(startnums1[i])+1,int(startnums1[i+1])+1):
             line=linecache.getline(filename6, j)
             store6.append(line)

   
   store7=[]
   for item2 in uniquescaffolds2:
      if (item in item2) and (item2 in item):
         i=uniquescaffolds2.index(item2)
         for j in range(int(startnums2[i])+1,int(startnums2[i+1])+1):
             line=linecache.getline(filename7, j)
             store7.append(line)
   store8=[]
   for item2 in uniquescaffolds3:
      if (item in item2) and (item2 in item):
         i=uniquescaffolds3.index(item2)
         for j in range(int(startnums3[i])+1,int(startnums3[i+1])+1):
             line=linecache.getline(filename8, j)
             store8.append(line)
   store=store6+store7+store8
   store=set(store)
   store=list(store)
   with open("temp_" + filename1,'w') as ijk:
    for item in store:
      ijk.write("%s" %(item))
   ijk.close()
   reader = csv.reader(open("temp_" + filename1), delimiter="\t")
   reader=list(reader)
   reader = sorted(reader, key=lambda abc: int(abc[3]))
   #reader.sorted(key=itemgetter(3))
   for row0,row1,row2,row3,row4,row5,row6,row7,row8 in reader:
     #print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(row0,row1,row2,row3,row4,row5,row6,row7,row8)
     xyz.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(str(row0),str(row1),str(row2),str(row3),str(row4),str(row5),str(row6),str(row7),str(row8)))
     
xyz.close()
         
with open("output_" + filename1,'w') as xyz:
  xyz.write("%s%i%s" %("The number of matches between de novo and repbase:",m,"\n"))
  xyz.write("%s%i%s%s%f%s"  %("Number of de novo hits:",len(reader1),"\n","Percentage of de novo matches:",(m*100)/len(reader1),"\n"))
  xyz.write("%s%i%s%s%f%s" %("Number of repbase hits:",len(reader2),"\n","Percentage of repbase matches:",(m*100)/len(reader2),"\n"))

     
xyz.close() 

# Deleting temporary files
os.remove("repbase_" + filename1)
os.remove("new_" + filename1 + ".align")
os.remove("known_" + filename1)
os.remove("unknown_" + filename1)      
os.remove("repbasedenovo_" + filename1)
os.remove("denovo_" + filename1)          
os.remove("denovo2_" + filename1)
os.remove("repbaseonly_" + filename1)
os.remove("all_" + filename1)
os.remove("temp_" + filename1)

time_end = time.time()
print "Time taken to run this parser is:",time_end-time_start," seconds"
