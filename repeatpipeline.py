from __future__ import division
from collections import Counter
from os import system
from scipy import stats
import operator
import csv
import re
import sys
from Bio import SeqIO  ## Try to change this
import os
import commands
import numpy as np
import argparse
import linecache
import random
import datetime
import time
from sklearn import svm
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score

def modharv(filename6,filename7,filename8):
  #filename6="element.gff" # Repeatmodeler gff file
  #filename7="harvest.gff" # LTR harvest gff file
  #filename8="Pila.fa" # Fasta file

  fastanames=[]
  fastalines=[]
  i=1
  for row in open(os.path.join("../repeatmasker", filename8), "r"):
   if row.startswith(">") == True:
       m=row.find(" ")
       if m!=-1:
         fastanames.append(row[1:].rstrip("\r\n")[:m-1])
       else:
         fastanames.append(row[1:].rstrip("\r\n"))
       #fastanames.append(row[1:])
       fastalines.append(i)
   i+=1
  fastalines.append(i)
  #print fastanames

  with open("elementtemp.gff",'w') as xyz:
    for line in open(os.path.join("../repeatmasker/full", filename6), "r"):
      if "##" not in line:
         xyz.write(line)
  xyz.close()
  with open("harvesttemp.gff",'w') as xyz:
    for line in open(filename7, "r"):
      if "###" not in line:
         xyz.write(line)
  xyz.close()
  with open("elementtemp.gff",'r') as f1:
        numberlines1=len(f1.readlines())
  f1.close()
  with open("harvesttemp.gff",'r') as f2:
        numberlines2=len(f2.readlines())
  f2.close()
  #print "numberlines1:",numberlines1,"numberlines2:",numberlines2
  elements=0
  genes=0
  startnums1=[]
  fastanames1=[]
  uniquescaffolds1=[]
  location1=[]
  row4=''
  m=-1
  strt=0
  for line1 in range(numberlines1):
       row1=linecache.getline("elementtemp.gff",line1+1)
       row1=row1.split("\t")
       m+=1       
       if len(row1)>5:
           if strt==0:
              temp=m
              strt=1  
           fastanames1.append(row1[0])
           try:
             if str(row1[0])!=str(row4[0]):
                location1.append(m)
                uniquescaffolds1.append(row1[0])
           except IndexError:
             location1.append(m)
             uniquescaffolds1.append(row1[0])
           row4=row1
  location1.append(numberlines1)
  #print "Fulls:",len(uniquescaffolds1),len(location1)
  sys.stdout.flush()
  startnums1=np.zeros(len(uniquescaffolds1)+2)
  startnums1[0]=temp
  i=1
  counts1=Counter(fastanames1)
  for item in uniquescaffolds1:
           startnums1[i]=counts1[item]+temp
           temp=startnums1[i]
           i+=1
  startnums1[len(uniquescaffolds1)]=len(fastanames1)
  #print "Uniquescaffolds1",uniquescaffolds1
  #print "Startnums1",startnums1
  sys.stdout.flush()

  startnums2=[]
  store7=[]
  fastanames2=[]
  uniquescaffolds2=[]
  location2=[]
  seq=[]
  fasta=[]
  row3=''
  m=-1
  strt=0
  totalhh=0
  for line2 in range(numberlines2):
       m+=1
       row2=linecache.getline("harvesttemp.gff",line2+1)  
       if (row2.startswith("###")==True) or (row2.startswith("##gff-version")==True): 
           continue
       elif row2.startswith("##sequence-region")==True:
           row2=row2.split(" ") 
           seq.append(row2[3])
       elif row2.startswith("#")==True:
           row2=row2.split(" ")
           fasta.append(row2[0][1:].rstrip("\r\n"))
       elif row2.startswith("seq")==True:
           if strt==0:
              temp=m
              strt=1
           row2=row2.split("\t") 
           fastanames2.append(row2[0])
           try:
             if str(row2[0])!=str(row3[0]):
               #print "row2:",row2,"\n","row3:",row3
               location2.append(m)

               uniquescaffolds2.append(row2[0])
           except IndexError:
             location2.append(m)
             uniquescaffolds2.append(row2[0])
           row3=row2

  location2.append(numberlines2)
  startnums2=np.zeros(len(uniquescaffolds2)+2)
  startnums2[0]=temp
  i=1
  #print "startnums2",len(startnums2)
  counts2=Counter(fastanames2)
  for item in uniquescaffolds2:
           startnums2[i]=counts2[item]+temp
           temp=startnums2[i]
           i+=1
  startnums2[len(uniquescaffolds2)]=len(fastanames2)
  #print "Uniquescaffolds2",uniquescaffolds2
  #print "Startnums2",startnums2
  #print "Seq",seq,"\n","Fasta",fasta

  with open(filename6[:len(filename6)-12] + "_linked.gff",'w') as ijk:
   with open(filename6[:len(filename6)-12] + "_alignments.gff",'w') as abc:
    with open(filename6[:len(filename6)-12] + "_linkedpartial.gff",'w') as pqr:
      for item2 in uniquescaffolds1:
       j=uniquescaffolds1.index(item2)
       for i in range(int(startnums1[j])+1,int(startnums1[j+1])+1):
         line=linecache.getline("elementtemp.gff", i)
         elements+=1  
         line3=line.split("\t")
         try:
              l=fasta.index(item2)
              m=seq[l]
         except ValueError:
              continue
         
         for k in range(int(startnums2[l])+1,int(startnums2[l+1])+1): 
           line2=linecache.getline("harvesttemp.gff", k)
           #print line,fasta[l],"\n",line2
           line4=line2.split("\t")
           if (len(line3) > 5) and (len(line4) > 5):
            if line4[2]=="repeat_region":
             if (int(line4[3]) < int(line3[3])) and (int(line3[4]) < int(line4[4])): # Gene
                if (int(line3[4])-int(line3[3]) >= 80) and (int(line3[4])-int(line3[3]) >= 0.8*(int(line4[4])-int(line4[3]))): # Invoking the 80/80 rule
                     ijk.write("%s%s" %(line,line2))
                     abc.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(fasta[l],line4[1],line4[2],line4[3],line4[4],line4[5],line4[6],line4[7],line4[8]))
                     sys.stdout.flush()
                     genes+=1
                else:
                     pqr.write("%s%s" %(line,line2))
             if (int(line3[3]) < int(line4[3])) and (int(line4[4]) < int(line3[4])): # Gene
                if (int(line4[4])-int(line4[3]) >= 80) and (int(line4[4])-int(line4[3]) >= 0.8*(int(line3[4])-int(line3[3]))): # Invoking the 80/80 rule
                     ijk.write("%s%s" %(line,line2))
                     abc.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(fasta[l],line4[1],line4[2],line4[3],line4[4],line4[5],line4[6],line4[7],line4[8]))
                     sys.stdout.flush()
                     genes+=1
                else:
                     pqr.write("%s%s" %(line,line2))                     
             #if (int(line4[3]) < int(line3[3])) and (int(line4[4]) < int(line3[4])): # Gene
              #pqr.write("%s%s" %(line,line2))
              #sys.stdout.flush()
              #genes+=1          
             #if (int(line3[3]) < int(line4[3])) and (int(line3[4]) < int(line4[4])): # Gene
              #pqr.write("%s%s" %(line,line2))
              #sys.stdout.flush()
              #genes+=1       
  pqr.close()
  abc.close()
  ijk.close()

  with open(filename6[:len(filename6)-12] + "_linked.gff",'r') as f3:
        numberlines3=len(f3.readlines())
  f3.close()
  rmod=[]
  family=[]
  harvest=[]

  for i in range(numberlines3):
      line=linecache.getline(filename6[:len(filename6)-12] + "_linked.gff",i+1)
      line=line.split()
      if line[1]=="RepeatMasker":
         line2=linecache.getline(filename6[:len(filename6)-12] + "_linked.gff",i+2)
         line2=line2.split()
         a=line[8].find("element=")
         harv=line[8][a+8:]
         b=harv.find(";")
         c=line[8].find("family=")
         harv2=line[8][c+7:]
         d=harv2.find(";")
         rmod.append(harv[:b])
         family.append(harv2[:d])
         harvest.append(line2[0])
  uniquermod=set(rmod)
  uniquermod=list(uniquermod)

  with open("modharv.txt",'w') as xyz:
   xyz.write("Rmod\tFamily\tLTR Harvest\n")
   for mod in uniquermod:
      hits=[]
      i=-1      
      for mod2 in rmod:
         i+=1
         if str(mod)==str(mod2):
            hits.append(harvest[i])
            fam=family[i]
      xyz.write("%s\t%s\t" %(mod,fam))
      xyz.write(' '.join(hits))
      xyz.write('\n')
  xyz.close()

  # Removing files
  os.remove("elementtemp.gff")
  os.remove("harvesttemp.gff")

  print "No. of lines analyzed:",i,"elements:",elements,"Genes:",genes,"\n","No. of repeats identified",len(uniquermod)  


def comparegff(filename1,minlen,minid,reppath,repfile,classfolder,classified):
  filename5=open(os.path.join("../", filename1), "r")
  filename4=(filename1[:len(filename1)-8] + ".align")   
  filename12=open(os.path.join("../", filename4), "r")

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
      if alignments[i]>=minlen:
         highestidentity[i]=percentidentity[i]

  # Extracting repeats which obey the 80/80 rule

  i=-1
  unknowns=0
  with open("known_" + filename1,'w') as xyz:
   with open("unknown_" + filename1,'w') as pqr:
    for line in open("all_" + filename1,'r'):
      i+=1
      if (alignments[i] >= minlen) and (percentidentity[i] >= minid):  # Invoking the 80/80 rule
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
  filename15="../../repeatmodeler/" + classfolder + "/" + classified


  items=[]
  conns=[]
  for seq_record in SeqIO.parse(filename15, "fasta"):
       identity=(seq_record.id)
       hashh=identity.find("#")
       imp=identity[hashh+1:len(identity)]
       conn=identity[:hashh]
       conns.append(conn)
       items.append(imp)

  lenfile=0
  for line in open("repbase_" + filename1,'r'):
     line=line.split()
     if len(line)>0:
       lenfile+=1
  reader2=np.loadtxt("repbase_" + filename1, delimiter='\t',dtype='str')
  #print "Length of reader2:", len(reader2)
  i=0
  detectors2=[]
  if lenfile==1:
     interval=np.zeros(1)
     for line in open("repbase_" + filename1,'r'):
         line=line.split("\t")
         s=line[8]
         s2=s[14:]
         s3=s2.find('"')
         detectors2.append(s2[:s3])
         i+=1       
  else:
     interval=np.zeros(len(reader2))
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

  filename21=reppath + repfile
  filename22=open(os.path.join(reppath, repfile), "r")

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
  if minlen==0 and minid==0:
     namefile=filename1[:len(filename1)-8] + "_final0000.gff"
  else:
     namefile=filename1[:len(filename1)-8] + "_final8080.gff"
  with open(namefile,'w') as xyz:
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
   try:  
    xyz.write("%s%i%s%s%f%s"  %("Number of de novo hits:",len(reader1),"\n","Percentage of de novo matches:",(m*100)/len(reader1),"\n"))
    xyz.write("%s%i%s%s%f%s" %("Number of repbase hits:",len(reader2),"\n","Percentage of repbase matches:",(m*100)/len(reader2),"\n"))
   except ZeroDivisionError:
    xyz.write("%s%i\n"  %("Number of de novo hits:",len(reader1)))
    xyz.write("%s%i\n" %("Number of repbase hits:",len(reader2)))   
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

def elementext(filename6,filename7,speciesname):
  #filename6="Pilafosfull_final8080.gff"
  #filename7="Pila.fa"
  reader1=np.loadtxt(filename6, delimiter='\t',dtype='str')
  fastanames=[]
  fastalines=[]
  i=1
  for row in open(os.path.join("../", filename7), "r"):
   if row.startswith(">") == True:
       m=row.find(" ")
       if m!=-1:
         fastanames.append(row[1:].rstrip("\r\n")[:m-1])
       else:
         fastanames.append(row[1:].rstrip("\r\n")) 
       fastalines.append(i)
   i+=1
  fastalines.append(i)
  seqlen=[]
  seq=[]
  #seq=[[] for i in range(len(fastanames))]
  s=-1
  #print len(seq),len(seqlen),len(fastanames)
  for name in fastanames:
     s+=1
     seq.append(name)
     k=fastanames.index(name)
     entry1=fastalines[k]
     entry2=fastalines[k+1]
     sums=0
     for m in range(entry1+1,entry2):
       abc=linecache.getline("../"+filename7,m)
       sums+=len(abc)
     seqlen.append(sums)
  k=0
  prev=""
  #print fastanames
  with open(filename6[:len(filename6)-14] + "_element.gff",'w') as xyz:
   xyz.write("##gff-version 3\n")
   for i in range(len(reader1)):
         line=linecache.getline(filename6,i+1)
         row0,row1,row2,row3,row4,row5,row6,row7,row8=line.split("\t")
         k+=1
         if row0 != prev:
            num=fastanames.index(row0)
            end=seqlen[num]
            xyz.write("%s   %s %s %s\n" %("##sequence-region",row0,"1",end))
            prev=row0
         column=row8
         element=column[:column.find("(")-1]
         family=column[column.find("(")+1:column.find(")")]
         if family.find("-") != -1:
            common=family[family.find("-")+1:]
            family=family[:family.find("-")]
         else:
            common=""     
         if column.find("[") != -1:
            species=column[column.find("[")+1:column.find("]")]
         else:
             if speciesname != None: 
                species=speciesname
         if family=="Gypsy":
             family="LTR/Gypsy"
         if family=="Copia":
             family="LTR/Copia"                   
         if common=="":
           if "Simple_repeat" in row8:
              pass
           elif "Satellite" in row8:
              pass
           else:     
              xyz.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%i%s%s%s%s%s%s\n" %(str(row0),str(row1),"SO:0000726",str(row3),str(row4),str(row5),str(row6),str(row7),"ID=",k,";element=",element,";family=",family,";species=",species))  
         else:
           if "Simple_repeat" in row8:
              pass
           elif "Satellite" in row8:
              pass
           else:     
              xyz.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%i%s%s%s%s%s%s%s%s%s\n" %(str(row0),str(row1),"SO:0000726",str(row3),str(row4),str(row5),str(row6),str(row7),"ID=",k,";element=",common,"(",element,")",";family=",family,";species=",species))  

time_start = time.time()
parser = argparse.ArgumentParser(
     prog='repeatpipeline',
     usage='''Usage: python repeatpipeline.py {options}
Necessary:
--fasta [Name of fasta file]
Optional: 
--name [name of genome (optional)] 
--path [Path to fasta file (optional)] 
--cpus [Number of CPU's (optional)] 
--pog [Percentage of genome for de novo repeat detection (optional)] 
--minlength [Minimum length to be considered as full-length] 
--minidentity [Minimum percent identity (0 to 100) to be considered as full-length]
--repbase [Reference database (from repbase) (optional)]
--repbasepath [Path to reference database, by default it assumes the location of this script  (optional)]
--LTRharvest [Type 0 to skip running LTR harvest]
--LTRharvestpath [Directory of LTR harvest path (optional)]''',
     description='''This program is a complete pipeline for identifying and characterizing transposable elements in a given genome.''',
     epilog='''It requires numpy, matplotlib, scikit and biopython libraries for execution.''')
parser.add_argument('--fasta', type=str, help='The query fasta file', required=True)
parser.add_argument('--name', type=str, help='The name of genome (optional)', required=False)
parser.add_argument('--pog', type=int, help='Percentage of genome for de novo repeat detection, please give whole numbers (optional)', required=False)
parser.add_argument('--path', type=str, help='Path to input files [Not needed if genome file in same directory as this script]', required=False)
parser.add_argument('--cpus', type=int, help="Number of CPU's (optional)", required=False)
parser.add_argument('--LTRharvest', type=int, help="Type 0 to skip running LTR harvest", required=False)
parser.add_argument('--minlength', type=int, help="Minimum length to be considered as full-length (default 80bp)", required=False)
parser.add_argument('--minidentity', type=float, help="Minimum percent identity (0 to 100) to be considered as full-length (default 80%)", required=False)
parser.add_argument('--LTRharvestpath', type=str, help="Directory of LTR harvest path (optional)", required=False)
parser.add_argument('--repbase', type=str, help='Reference database (from repbase) (optional)', required=False)
parser.add_argument('--repbasepath', type=str, help='Path to reference database, by default it assumes the location of this script  (optional)', required=False)
args = parser.parse_args()

filename1=args.fasta
path=args.path
nameofg=args.name
cpus=args.cpus
pog=args.pog
ltrharveststatus=args.LTRharvest
ltrharvestpath=args.LTRharvestpath
repbasefile=args.repbase
minlength=args.minlength
minidentity=args.minidentity
repbasepath=args.repbasepath
if path == None:
   filename5=open(filename1,"r")
else:
   filename5=open(os.path.join(path, filename1), "r")

if path == None:
      filename6=filename1
else:
      if path[len(path)-1]=="/":
                filename6=path + filename1
      else:
                filename6=path + "/" + filename1

if nameofg == None:
   nameofg="temp"

if cpus == None:
   cpus=1    

# Making a new directory
status, homedir = commands.getstatusoutput("pwd")
if repbasepath == None:
    repbasepath=homedir
if repbasepath[len(repbasepath)-1]=="/":
   pass
else:
   repbasepath=repbasepath + "/"
foldername="repeat_" + nameofg + "_" + str(datetime.date.today()) + "_"  + str(time.time())
print "Name of folder:",foldername
os.system("mkdir " + foldername)
os.chdir(foldername)

# Running repeatmodeler
os.system("mkdir repeatmodeler")
os.chdir("repeatmodeler")
if pog == None:
   if path == None:
        os.system("BuildDatabase -name " + nameofg + " -engine ncbi " + filename6)
   else:
        os.system("BuildDatabase -name " + nameofg + " -engine ncbi " + filename6)
   os.system("RepeatModeler -database " + nameofg + " -engine ncbi -pa " + str(cpus))
   repeatmodelertime = time.time()
   print "RepeatModeler run successful.\n Time required to run repeatmodeler is:",repeatmodelertime-time_start

else:
   if path==None:
      filename8=open(filename1,'r')
   else:
      filename8=open(os.path.join(path, filename1), "r")

   if path == None:
      filename9=filename1
   else:
      if path[len(path)-1]=="/": 
         filename9=path + filename1
      else:
         filename9=path + "/" + filename1

   fastanames1=[]
   fastalines1=[]
   i=1
   for row in filename8:
     if row.startswith(">") == True:
        fastanames1.append(row[1:])
        fastalines1.append(i)
     i+=1
   fastalines1.append(i)
   # Calulate total number of bases

   genome=""

   for name in fastanames1:
       k=fastanames1.index(name)
       entry1=fastalines1[k]
       entry2=fastalines1[k+1]
       for m in range(entry1+1,entry2):
           genome+=linecache.getline(filename9,m)
   pieces=100 # Change this number to change the no. of pieces of the genome 
   splice=len(genome)/pieces  
   j=0
   splitnames=[]
   for i in range(pieces):
      c = filename1 + str(i+1) + ".fa"
      splitnames.append(c)
      breaks=0
      with open(c,'w') as ijk:        
         while (breaks < splice) and j < len(fastanames1):
             entry1=fastalines1[j]
             entry2=fastalines1[j+1]
             print entry1,entry2
             sequence=""
             for m in range(entry1+1,entry2):
                 sequence+=linecache.getline(filename9,m)
             ijk.write("%s%s\n%s\n" % (">",fastanames1[j].rstrip('\r\n'), sequence.rstrip('\r\n'))) # Making the fasta file
             #print j,fastanames1[j]
             breaks+=len(sequence)
             j+=1
      ijk.close()
   random.shuffle(splitnames)   
   tempset=""
   for i in range(pog):
      tempset+=splitnames[i]
      tempset+=" "
   os.system("cat " + tempset + "> repeatmodelertrainingset.fasta")
   print "Name of database:", nameofg,"\n"
   os.system("BuildDatabase -name " + nameofg + " -engine ncbi repeatmodelertrainingset.fasta")
   print "Database built"
   os.system("RepeatModeler -database " + nameofg + " -engine ncbi -pa " + str(cpus))
   repeatmodelertime = time.time()
   print "RepeatModeler run successful.\n Time required to run repeatmodeler is:",repeatmodelertime-time_start

status, output = commands.getstatusoutput("ls")

output2=output.split(" ")
output3=output.split("\n")
for word in output3:
   if word.startswith("RM_") == True:
      resultsfolder=word

os.chdir(resultsfolder)

######## Repeat classification starts

# Running CENSOR

if repbasefile == None:
    os.system("cp " + repbasepath + "all_plantupdated.ref all_plantupdated.ref")
    repbasefile="all_plantupdated.ref"
else:
    os.system("cp " + repbasepath + repbasefile + " " + repbasefile)
os.system("censor.ncbi consensi.fa.classified -lib " + repbasefile  + " -tab -no_simple -s")

# Incorporating CENSOR results into the main repeatmodeler library

with open("unknown8080.txt",'w') as xyz:
 for line in open("consensi.fa.classified.map",'r'):
    if "Unknown" in line:
       line2=line.split()
       if float(line2[7]) > 0.8 and (int(line2[2])-int(line2[1]) >= 80): # Applying the 80/80 cut-off.
          xyz.write("%s" %(line))
xyz.close()

# Incorporating all de novo results into a combined de novo library.

with open(nameofg + "_finalrepeatclassifications.fasta",'w') as xyz:
 for line2 in open("consensi.fa.classified",'r'):
   for line in open("unknown8080.txt",'r'):
        m=line.find("#")
        n=line2.find("#")
        if m != -1:
           phrase1=line[:m]
        if n != -1:
           phrase2=line2[1:n]
        if (phrase1 in phrase2) and (phrase2 in phrase1) and (m != -1) and (n != -1):
           line=line.split()
           answer=line[3]
           if ("Gypsy" in answer) or ("GYPSY" in answer) or ("gypsy" in answer) or ("PtOuachita" in answer) or ("PtBastrop" in answer) or ("PtOzark" in answer) or ("PtAppalachian" in answer) or ("PtAngelina" in answer) or ("PtTalladega" in answer)  or ("IFG7" in answer): 
              answer2 ="LTR/Gypsy-"+answer
           elif ("Copia" in answer) or ("COPIA" in answer) or ("copia" in answer) or ("PtCumberland" in answer) or ("PtPineywoods" in answer) or ("PtConagree" in answer) or ("Silava" in answer):
              answer2="LTR/Copia-"+answer
           elif ("LINE" in answer) :
              answer2="LINE/" + answer
           elif ("SINE" in answer) :
              answer2="SINE/"+answer
           elif ("Helitron" in answer) or ("helitron" in answer):
              answer2="DNA/Helitron-"+answer
           elif ("Harbinger" in answer) or ("harbinger" in answer):
              answer2="DNA/Harbinger-"+answer
           elif ("hAT" in answer):
              answer2="DNA/hAT" +answer
           elif ("MuDR" in answer):
              answer2="DNA/MULE-"+answer
           elif ("EnSpm" in answer):
              answer2="DNA/CMC-"+answer
           elif ("DNA" in answer) :
              answer2="DNA/"+answer
           else:
              answer2=answer
           answer2=answer2.strip("\r\n")
           o=line2.find("(")
           q=line2.find(")")
           recon=line2[o:q+1]
           print "CENSOR",answer2   
           line2=phrase1 + "#" + str(answer2) + " " + recon + " CENSOR" + "\n"
   xyz.write("%s" %(line2))        
xyz.close()
os.system("cat " + nameofg + "_finalrepeatclassifications.fasta " + repbasefile + " > " + nameofg + "_completelib.fasta")
status, repeatmodelerdir = commands.getstatusoutput("pwd")
os.chdir("..")
os.chdir("..")

# Making 10MB size bins for repeatmasker and LTR harvest
os.system("mkdir repeatmasker LTRharvest")
os.chdir("repeatmasker")
status, repeatmaskerdir = commands.getstatusoutput("pwd")
repeatmaskerstart=time.time()

genomefilesize=os.path.getsize(filename6) # Checking file-size of genome
eachbin=10485760 # 10MB each bin
numbins=int(genomefilesize)//eachbin
numbins+=1

splice=len(genome)/numbins  
j=0
splitnames=[]
for i in range(pieces):
      c = filename1 + str(i+1) + ".fa"
      splitnames.append(c)
      breaks=0
      with open(c,'w') as ijk:        
         while (breaks < splice) and j < len(fastanames1):
             entry1=fastalines1[j]
             entry2=fastalines1[j+1]
             print entry1,entry2
             sequence=""
             for m in range(entry1+1,entry2):
                 sequence+=linecache.getline(filename9,m)
             ijk.write("%s%s\n%s\n" % (">",fastanames1[j].rstrip('\r\n'), sequence.rstrip('\r\n'))) # Making the fasta file
             #print j,fastanames1[j]
             breaks+=len(sequence)
             j+=1
      ijk.close()

# Removing blank files

for i in range(pieces):
   c = filename1 + str(i+1) + ".fa"
   if int(os.path.getsize(c))==0:
      os.system("rm " + c)
      pieces-=1

print "Total number of bins:",pieces

# Running Repeatmasker and LTR harvest

for i in range(pieces):
    os.system("RepeatMasker " + str(filename1) + str(i+1) + ".fa"  + " -dir " + str(repeatmaskerdir) + " -lib " + str(repeatmodelerdir) + "/" + str(nameofg) + "_completelib.fasta" + " -pa " + str(cpus) + " -gff -a -noisy -low")
    if (ltrharveststatus==None) or (ltrharveststatus==1):
      ltrharveststart=time.time()
      os.chdir("../LTRharvest")
      if ltrharvestpath==None:
         os.system("gt suffixerator -db " + "../repeatmasker/" + str(filename1) + str(i+1) + ".fa" + " -suf -lcp")
         os.system("gt -j " + str(cpus) + " ltrharvest -index " + str(filename1) + str(i+1) + ".fa" + " -gff3 " + "harvest" + str(i+1) + ".gff")
      else:
         os.system(ltrharvestpath + "/" + "gt suffixerator -db " + "../repeatmasker/" + str(filename1) + str(i+1) + ".fa" + " -suf -lcp")
         os.system(ltrharvestpath + "/" + "gt -j " + str(cpus) + " ltrharvest -index " + str(filename1) + str(i+1) + ".fa" + " -gff3 " + "harvest" + str(i+1) + ".gff")
      os.chdir("../repeatmasker")   
    elif ltrharveststatus==0:
      os.mknod("harvest" + str(i+1) + ".gff")

print "Time taken to run RepeatMasker and LTR harvest is: ",time.time() - repeatmaskerstart      

# Parsing out full-length and partial-length repeats

os.system("mkdir full")
os.system("mkdir partial")
os.chdir("full")
status, fulllengthdir = commands.getstatusoutput("pwd")

if minlength==None:
   minlength=80
if minidentity==None:
   minidentity=80

for i in range(pieces):
    comparegff(filename1 + str(i+1) + ".fa.out.gff",minlength,minidentity,repbasepath,repbasefile,resultsfolder,nameofg + "_finalrepeatclassifications.fasta") # Full-length repeats
    elementext(filename1 + str(i+1) + ".fa_final8080.gff",filename1 + str(i+1) + ".fa",nameofg)
    os.chdir("../../LTRharvest")
    modharv(filename1 + str(i+1) + ".fa_element.gff","harvest" + str(i+1) + ".gff",filename1 + str(i+1) + ".fa")
    os.chdir("../repeatmasker/partial")   
    comparegff(filename1 + str(i+1) + ".fa.out.gff",0,0,repbasepath,repbasefile,resultsfolder,nameofg + "_finalrepeatclassifications.fasta") # (Full + partial)-length repeats
    elementext(filename1 + str(i+1) + ".fa_final0000.gff",filename1 + str(i+1) + ".fa",nameofg)
    os.chdir("../full")
