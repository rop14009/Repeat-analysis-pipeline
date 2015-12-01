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
         fastanames.append(row[1:].replace("\n","")[:m-1])
       else:
         fastanames.append(row[1:].replace("\n",""))
       #fastanames.append(row[1:])
       fastalines.append(i)
   i+=1
  fastalines.append(i)
  #print fastanames

  with open("element"+filename6+"temp.gff",'w') as xyz:
    for line in open(os.path.join("../repeatmasker/full", filename6), "r"):
      if "##" not in line:
         xyz.write(line)
  xyz.close()
  with open("harvest"+filename7+"temp.gff",'w') as xyz:
    for line in open(filename7, "r"):
      if "###" not in line:
         xyz.write(line)
  xyz.close()
  with open("element"+filename6+"temp.gff",'r') as f1:
        numberlines1=len(f1.readlines())
  f1.close()
  with open("harvest"+filename7+"temp.gff",'r') as f2:
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
       row1=linecache.getline("element"+filename6+"temp.gff",line1+1)
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
       row2=linecache.getline("harvest"+filename7+"temp.gff",line2+1)  
       if (row2.startswith("###")==True) or (row2.startswith("##gff-version")==True): 
           continue
       elif row2.startswith("##sequence-region")==True:
           row2=row2.split(" ") 
           seq.append(row2[3])
       elif row2.startswith("#")==True:
           row2=row2.split(" ")
           fasta.append(row2[0][1:].replace("\n",""))
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
         line=linecache.getline("element"+filename6+"temp.gff", i)
         elements+=1  
         line3=line.split("\t")
         try:
              l=fasta.index(item2)
              m=seq[l]
         except ValueError:
              continue
         try:
          for k in range(int(startnums2[l])+1,int(startnums2[l+1])+1): 
           line2=linecache.getline("harvest"+filename7+"temp.gff", k)
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
         except IndexError:
              print "l:",l,'\n',"Length of startnums2:",len(startnums2),"\n","startnums2:",startnums2
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
  os.remove("element"+filename6+"temp.gff")
  os.remove("harvest"+filename7+"temp.gff")

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
         fastanames.append(row[1:].replace("\n","")[:m-1])
       else:
         fastanames.append(row[1:].replace("\n","")) 
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
     usage='''python repeatpipeline.py {options}
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
parser.add_argument('--fasta', type=str, help='The genome fasta file', required=True)
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
#filename5.close()
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
   fastanames1=[]
   fastalines1=[]
   i=1
   for row in filename5:
     if row.startswith(">") == True:
        fastanames1.append(row[1:])
        fastalines1.append(i)
     i+=1
   fastalines1.append(i)
   # Calulate total number of bases
   filename5.close()
   genome=""

   for name in fastanames1:
       k=fastanames1.index(name)
       entry1=fastalines1[k]
       entry2=fastalines1[k+1]
       for m in range(entry1+1,entry2):
           genome+=linecache.getline(filename6,m)
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
                 sequence+=linecache.getline(filename6,m)
             ijk.write("%s%s\n%s\n" % (">",fastanames1[j].replace("\n",""), sequence.replace("\n",""))) # Making the fasta file
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

fastanames1=[]
fastalines1=[]
i=1

if path == None:
   filename5=open(filename1,"r")
else:
   filename5=open(os.path.join(path, filename1), "r")

for row in filename5:
  if row.startswith(">") == True:
     fastanames1.append(row[1:])
     fastalines1.append(i)
  i+=1
fastalines1.append(i)
filename5.close()

genome=""

for name in fastanames1:
    k=fastanames1.index(name)
    entry1=fastalines1[k]
    entry2=fastalines1[k+1]
    for m in range(entry1+1,entry2):
        genome+=linecache.getline(filename6,m)

splice=len(genome)/numbins

# Determining number of bins of genome required

if genomefilesize<10485760:  
   pieces=1
elif (genomefilesize>10485760) and (genomefilesize<=104857600):
   pieces=10
elif (genomefilesize>104857600) and (genomefilesize<=1048576000):   
   pieces=100
elif (genomefilesize>1048576000) and (genomefilesize<=10485760000):
   pieces=1000
elif (genomefilesize>10485760000):
   pieces=10000
   
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
                 sequence+=linecache.getline(filename6,m)
             ijk.write("%s%s\n%s\n" % (">",fastanames1[j].replace("\n",""), sequence.replace("\n",""))) # Making the fasta file
             #print j,fastanames1[j]
             breaks+=len(sequence)
             j+=1
      ijk.close()

# Removing blank files

for i in range(pieces):
   c = filename1 + str(i+1) + ".fa"
   if int(os.path.getsize(c))==0:
      os.remove(c)
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
    if ltrharveststatus==1:
       os.chdir("../../LTRharvest")
       modharv(filename1 + str(i+1) + ".fa_element.gff","harvest" + str(i+1) + ".gff",filename1 + str(i+1) + ".fa")
       os.chdir("../repeatmasker/partial")
    else:
       os.chdir("../partial")
    comparegff(filename1 + str(i+1) + ".fa.out.gff",0,0,repbasepath,repbasefile,resultsfolder,nameofg + "_finalrepeatclassifications.fasta") # (Full + partial)-length repeats
    elementext(filename1 + str(i+1) + ".fa_final0000.gff",filename1 + str(i+1) + ".fa",nameofg)
    os.chdir("../full")
os.system("cat *.fa_final8080.gff > " + nameofg + "full_final8080.gff")

os.chdir("../../LTRharvest")    
# Looking for unknown LTR elements detected by LTR harvest
unknownelements=[]
if ltrharveststatus==1:
   os.system("cat *_linked.gff > all_linked.gff")
   with open("LTRharvesthits.txt",'w') as xyz: 
     for line in open("all_linked.gff",'r'):
       if "RepeatMasker" in line:
          if "Unknown" in line:
             line2=line.split('\t')
             m=line2[8].find(";element=")
             temp1=line2[8][m+9:]
             n=line2[8].find(";family=")
             element=temp1[:n]
             xyz.write("%s\n" %(element))
             unknownelements.append(element)
   xyz.close()
   if len(unknownelements)>0:
    with open(nameofg + "_finalrepeatclassifications_LTRharvest.fasta",'w') as pqr:
     for row in open(os.path.join(repeatmodelerdir, nameofg + "_finalrepeatclassifications.fasta"), "r"):
      if row.startswith(">")==True:
        for element in unknownelements:  
                 hits1=re.findall('\\b' + str(element) + '\\b', row)
                 if len(hits1)>0:
                   row=row.replace("Unknown","LTR")
                   break
        pqr.write(row)  
      else:
        pqr.write(row)
    pqr.close()      
   else:
    os.system("cp " + repeatmodelerdir + "/" + nameofg + "_finalrepeatclassifications.fasta " + nameofg + "_finalrepeatclassifications_LTRharvest.fasta")
else:
 os.system("cp " + repeatmodelerdir + "/" + nameofg + "_finalrepeatclassifications.fasta " + nameofg + "_finalrepeatclassifications_LTRharvest.fasta")

# Its repbaseext2.py code 
 
#filename1=args.gff
#filename2=args.repbase
#filename5=args.denovo
#path=args.path
#name2=args.name

filename3=open(os.path.join("../repeatmasker/full", nameofg + "full_final8080.gff"), "r")
filename4=open(os.path.join(repbasepath, repbasefile), "r")
filename6=open(nameofg + "_finalrepeatclassifications_LTRharvest.fasta","r")
filename9=repbasepath + repbasefile
elements=[]
for line in filename3:
   line2=line.split()
   abc=line2[8]
   m=line.find('(')
   n=line.find(')')
   elements.append(abc[:m-1])
filename3.close()
uniqueelements=set(elements)
uniqueelements=list(uniqueelements)
uniqueelements=sorted(uniqueelements)  
   
numberlines=0
fastanames=[]
fastalines=[]
i=1
for row in filename4: # Repbase
 if row.startswith(">") == True:
     row2=row[1:]
     name=row2.split()
     fastanames.append(name[0])
     fastalines.append(i)
 i+=1   
fastalines.append(i)
filename4.close()
print "Repbase",fastanames[:10]
fastanames2=[]
fastalines2=[]
i=1
for row in filename6: # Denovo
 if row.startswith(">") == True:
     row2=row[1:]
     name=row2.split()
     fastanames2.append(name[0])
     fastalines2.append(i)
 i+=1   
fastalines2.append(i)
filename6.close()
print "Denovo",fastanames2[:10]
i=-1
total=0
counts=0
counts2=0
counts3=0
with open("library_" + nameofg + ".fasta",'w') as pqr:
 for name in fastanames2:
   n=fastanames2.index(name)
   entry1=fastalines2[n]
   entry2=fastalines2[n+1]
   seq1=linecache.getline(nameofg + "_finalrepeatclassifications_LTRharvest.fasta",entry1)
   a=seq1.find('#')
   b=seq1.find('(')
   elementname=seq1[1:a]
   familyname=seq1[a+1:b-1]
   if "Pt" in elementname:
      seq2=">"+familyname+"("+elementname+")"+"[Pinus taeda]"
   else:
      seq2=">"+familyname+"("+elementname+")"
   seq2=seq2.replace("/","-")
   pqr.write("%s\n" %(seq2))
   for m in range(entry1+1,entry2):
      pqr.write("%s" %(linecache.getline(nameofg + "_finalrepeatclassifications_LTRharvest.fasta",m)))

 for name in uniqueelements:
  if "rnd" in name:
   continue
  else:
   total+=1
   try:
      k=fastanames.index(name)
      entry1=fastalines[k]
      entry2=fastalines[k+1]
      seq1=linecache.getline(filename9,entry1)[1:]
      seq1=seq1.split("\t")
      elementname=seq1[0]
      familyname=seq1[1]
      speciesname=seq1[2]
      if (familyname == "Gypsy") or (familyname == "gypsy"):
         familyname="LTR-Gypsy"
      if (familyname == "Copia") or (familyname == "copia"):
         familyname="LTR-Copia"
      if (familyname == "Harbinger") or (familyname == "harbinger"):
         familyname="DNA-Harbinger"
      if (familyname == "Penelope") or (familyname == "penelope"):
         familyname="DNA-Penelope"
      if (familyname == "hAT"): 
         familyname="DNA-hAT"      
      if (familyname == "EnSpm"): 
         familyname="DNA-Enspm"
      if (familyname == "MuDR"): 
         familyname="DNA-MuDR"            
      seq2=">"+familyname+"-"+elementname+"-"+"["+speciesname.replace("\n","")+"]\n"    
      seq2=seq2.replace("rnd",nameofg + "_rnd")
      pqr.write("%s" %(seq2))
      for m in range(entry1+1,entry2):
         pqr.write("%s" %(linecache.getline(filename9,m)))    
      counts3+=1
   except ValueError:   
      pass 

print "Total no. of repbase hits:",total,len(fastanames2)


# Machine learning starts

os.system("mkdir ../ML")
os.chdir("../ML")

# Differentiatating between DNA transposons and LTR retrotransposons

fastanames=[]
fastalines=[]
i=1
for row in open(os.path.join("../LTRharvest", "library_" + nameofg + ".fasta"), "r"):
 if row.startswith(">") == True:
     fastanames.append(row[1:])
     fastalines.append(i)
 i+=1
fastalines.append(i)

# Getting rid of simple repeats and rRNA.       

with open("itsp_" + "library_" + nameofg + ".fasta",'w') as xyz:
 for name in fastanames:
    if (("Satellite" in name) or ("rRNA" in name) or ("Simple_repeat" in name) or ("buffer" in name)  or ("snRNA" in name) or ("SAT" in name)):
        pass
    elif ("STOWAWAY" in name) or ("Sola" in name) or ("sola" in name) or ("Crypton" in name) or ("hAT" in name) or ("MuDR" in name) or ("HARB" in name) or ("EnSpm" in name) or ("Harbinger" in name) or ("harbinger" in name) or ("Helitron" in name) or ("helitron" in name) or ("DNA" in name) or ("R1" in name) or ("L1" in name) or ("LTR" in name) or ("Retro" in name) or ("SINE" in name) or ("PtCumberland" in name) or ("PtPineywoods" in name) or ("PtAppalachian" in name) or ("PtTalladega" in name) or ("PtAngelina" in name) or ("PtOzark" in name) or ("PtOuachita" in name) or ("PtBastrop" in name) or ("PtConagree" in name) or ("LINE" in name) or ("Copia" in name) or ("copia" in name) or ("COPIA" in name) or ("Gypsy" in name) or ("GYPSY" in name) or ("gypsy" in name): 
        k=fastanames.index(name)
        entry1=fastalines[k]
        entry2=fastalines[k+1]
        xyz.write("%s%s\n" %(">",name.replace("\n","")))
        for m in range(entry1+1,entry2):
           xyz.write("%s" %(linecache.getline("../LTRharvest/" + "library_" + nameofg + ".fasta",m).upper()))
xyz.close()

fastanames1=[]
fastalines1=[]
i=1
for row in open("itsp_" + "library_" + nameofg + ".fasta",'r'):
 if row.startswith(">") == True:
     fastanames1.append(row[1:])
     fastalines1.append(i)
 i+=1
fastalines1.append(i)

tetranucleotides =('AAAAA', 'AAAAT', 'AAAAC', 'AAAAG', 'AAATA', 'AAATT', 'AAATC', 'AAATG', 'AAACA', 'AAACT', 'AAACC', 'AAACG', 'AAAGA', 'AAAGT', 'AAAGC', 'AAAGG', 'AATAA', 'AATAT', 'AATAC', 'AATAG', 'AATTA', 'AATTT', 'AATTC', 'AATTG', 'AATCA', 'AATCT', 'AATCC', 'AATCG', 'AATGA', 'AATGT', 'AATGC', 'AATGG', 'AACAA', 'AACAT', 'AACAC', 'AACAG', 'AACTA', 'AACTT', 'AACTC', 'AACTG', 'AACCA', 'AACCT', 'AACCC', 'AACCG', 'AACGA', 'AACGT', 'AACGC', 'AACGG', 'AAGAA', 'AAGAT', 'AAGAC', 'AAGAG', 'AAGTA', 'AAGTT', 'AAGTC', 'AAGTG', 'AAGCA', 'AAGCT', 'AAGCC', 'AAGCG', 'AAGGA', 'AAGGT', 'AAGGC', 'AAGGG', 'ATAAA', 'ATAAT', 'ATAAC', 'ATAAG', 'ATATA', 'ATATT', 'ATATC', 'ATATG', 'ATACA', 'ATACT', 'ATACC', 'ATACG', 'ATAGA', 'ATAGT', 'ATAGC', 'ATAGG', 'ATTAA', 'ATTAT', 'ATTAC', 'ATTAG', 'ATTTA', 'ATTTT', 'ATTTC', 'ATTTG', 'ATTCA', 'ATTCT', 'ATTCC', 'ATTCG', 'ATTGA', 'ATTGT', 'ATTGC', 'ATTGG', 'ATCAA', 'ATCAT', 'ATCAC', 'ATCAG', 'ATCTA', 'ATCTT', 'ATCTC', 'ATCTG', 'ATCCA', 'ATCCT', 'ATCCC', 'ATCCG', 'ATCGA', 'ATCGT', 'ATCGC', 'ATCGG', 'ATGAA', 'ATGAT', 'ATGAC', 'ATGAG', 'ATGTA', 'ATGTT', 'ATGTC', 'ATGTG', 'ATGCA', 'ATGCT', 'ATGCC', 'ATGCG', 'ATGGA', 'ATGGT', 'ATGGC', 'ATGGG', 'ACAAA', 'ACAAT', 'ACAAC', 'ACAAG', 'ACATA', 'ACATT', 'ACATC', 'ACATG', 'ACACA', 'ACACT', 'ACACC', 'ACACG', 'ACAGA', 'ACAGT', 'ACAGC', 'ACAGG', 'ACTAA', 'ACTAT', 'ACTAC', 'ACTAG', 'ACTTA', 'ACTTT', 'ACTTC', 'ACTTG', 'ACTCA', 'ACTCT', 'ACTCC', 'ACTCG', 'ACTGA', 'ACTGT', 'ACTGC', 'ACTGG', 'ACCAA', 'ACCAT', 'ACCAC', 'ACCAG', 'ACCTA', 'ACCTT', 'ACCTC', 'ACCTG', 'ACCCA', 'ACCCT', 'ACCCC', 'ACCCG', 'ACCGA', 'ACCGT', 'ACCGC', 'ACCGG', 'ACGAA', 'ACGAT', 'ACGAC', 'ACGAG', 'ACGTA', 'ACGTT', 'ACGTC', 'ACGTG', 'ACGCA', 'ACGCT', 'ACGCC', 'ACGCG', 'ACGGA', 'ACGGT', 'ACGGC', 'ACGGG', 'AGAAA', 'AGAAT', 'AGAAC', 'AGAAG', 'AGATA', 'AGATT', 'AGATC', 'AGATG', 'AGACA', 'AGACT', 'AGACC', 'AGACG', 'AGAGA', 'AGAGT', 'AGAGC', 'AGAGG', 'AGTAA', 'AGTAT', 'AGTAC', 'AGTAG', 'AGTTA', 'AGTTT', 'AGTTC', 'AGTTG', 'AGTCA', 'AGTCT', 'AGTCC', 'AGTCG', 'AGTGA', 'AGTGT', 'AGTGC', 'AGTGG', 'AGCAA', 'AGCAT', 'AGCAC', 'AGCAG', 'AGCTA', 'AGCTT', 'AGCTC', 'AGCTG', 'AGCCA', 'AGCCT', 'AGCCC', 'AGCCG', 'AGCGA', 'AGCGT', 'AGCGC', 'AGCGG', 'AGGAA', 'AGGAT', 'AGGAC', 'AGGAG', 'AGGTA', 'AGGTT', 'AGGTC', 'AGGTG', 'AGGCA', 'AGGCT', 'AGGCC', 'AGGCG', 'AGGGA', 'AGGGT', 'AGGGC', 'AGGGG', 'TAAAA', 'TAAAT', 'TAAAC', 'TAAAG', 'TAATA', 'TAATT', 'TAATC', 'TAATG', 'TAACA', 'TAACT', 'TAACC', 'TAACG', 'TAAGA', 'TAAGT', 'TAAGC', 'TAAGG', 'TATAA', 'TATAT', 'TATAC', 'TATAG', 'TATTA', 'TATTT', 'TATTC', 'TATTG', 'TATCA', 'TATCT', 'TATCC', 'TATCG', 'TATGA', 'TATGT', 'TATGC', 'TATGG', 'TACAA', 'TACAT', 'TACAC', 'TACAG', 'TACTA', 'TACTT', 'TACTC', 'TACTG', 'TACCA', 'TACCT', 'TACCC', 'TACCG', 'TACGA', 'TACGT', 'TACGC', 'TACGG', 'TAGAA', 'TAGAT', 'TAGAC', 'TAGAG', 'TAGTA', 'TAGTT', 'TAGTC', 'TAGTG', 'TAGCA', 'TAGCT', 'TAGCC', 'TAGCG', 'TAGGA', 'TAGGT', 'TAGGC', 'TAGGG', 'TTAAA', 'TTAAT', 'TTAAC', 'TTAAG', 'TTATA', 'TTATT', 'TTATC', 'TTATG', 'TTACA', 'TTACT', 'TTACC', 'TTACG', 'TTAGA', 'TTAGT', 'TTAGC', 'TTAGG', 'TTTAA', 'TTTAT', 'TTTAC', 'TTTAG', 'TTTTA', 'TTTTT', 'TTTTC', 'TTTTG', 'TTTCA', 'TTTCT', 'TTTCC', 'TTTCG', 'TTTGA', 'TTTGT', 'TTTGC', 'TTTGG', 'TTCAA', 'TTCAT', 'TTCAC', 'TTCAG', 'TTCTA', 'TTCTT', 'TTCTC', 'TTCTG', 'TTCCA', 'TTCCT', 'TTCCC', 'TTCCG', 'TTCGA', 'TTCGT', 'TTCGC', 'TTCGG', 'TTGAA', 'TTGAT', 'TTGAC', 'TTGAG', 'TTGTA', 'TTGTT', 'TTGTC', 'TTGTG', 'TTGCA', 'TTGCT', 'TTGCC', 'TTGCG', 'TTGGA', 'TTGGT', 'TTGGC', 'TTGGG', 'TCAAA', 'TCAAT', 'TCAAC', 'TCAAG', 'TCATA', 'TCATT', 'TCATC', 'TCATG', 'TCACA', 'TCACT', 'TCACC', 'TCACG', 'TCAGA', 'TCAGT', 'TCAGC', 'TCAGG', 'TCTAA', 'TCTAT', 'TCTAC', 'TCTAG', 'TCTTA', 'TCTTT', 'TCTTC', 'TCTTG', 'TCTCA', 'TCTCT', 'TCTCC', 'TCTCG', 'TCTGA', 'TCTGT', 'TCTGC', 'TCTGG', 'TCCAA', 'TCCAT', 'TCCAC', 'TCCAG', 'TCCTA', 'TCCTT', 'TCCTC', 'TCCTG', 'TCCCA', 'TCCCT', 'TCCCC', 'TCCCG', 'TCCGA', 'TCCGT', 'TCCGC', 'TCCGG', 'TCGAA', 'TCGAT', 'TCGAC', 'TCGAG', 'TCGTA', 'TCGTT', 'TCGTC', 'TCGTG', 'TCGCA', 'TCGCT', 'TCGCC', 'TCGCG', 'TCGGA', 'TCGGT', 'TCGGC', 'TCGGG', 'TGAAA', 'TGAAT', 'TGAAC', 'TGAAG', 'TGATA', 'TGATT', 'TGATC', 'TGATG', 'TGACA', 'TGACT', 'TGACC', 'TGACG', 'TGAGA', 'TGAGT', 'TGAGC', 'TGAGG', 'TGTAA', 'TGTAT', 'TGTAC', 'TGTAG', 'TGTTA', 'TGTTT', 'TGTTC', 'TGTTG', 'TGTCA', 'TGTCT', 'TGTCC', 'TGTCG', 'TGTGA', 'TGTGT', 'TGTGC', 'TGTGG', 'TGCAA', 'TGCAT', 'TGCAC', 'TGCAG', 'TGCTA', 'TGCTT', 'TGCTC', 'TGCTG', 'TGCCA', 'TGCCT', 'TGCCC', 'TGCCG', 'TGCGA', 'TGCGT', 'TGCGC', 'TGCGG', 'TGGAA', 'TGGAT', 'TGGAC', 'TGGAG', 'TGGTA', 'TGGTT', 'TGGTC', 'TGGTG', 'TGGCA', 'TGGCT', 'TGGCC', 'TGGCG', 'TGGGA', 'TGGGT', 'TGGGC', 'TGGGG', 'CAAAA', 'CAAAT', 'CAAAC', 'CAAAG', 'CAATA', 'CAATT', 'CAATC', 'CAATG', 'CAACA', 'CAACT', 'CAACC', 'CAACG', 'CAAGA', 'CAAGT', 'CAAGC', 'CAAGG', 'CATAA', 'CATAT', 'CATAC', 'CATAG', 'CATTA', 'CATTT', 'CATTC', 'CATTG', 'CATCA', 'CATCT', 'CATCC', 'CATCG', 'CATGA', 'CATGT', 'CATGC', 'CATGG', 'CACAA', 'CACAT', 'CACAC', 'CACAG', 'CACTA', 'CACTT', 'CACTC', 'CACTG', 'CACCA', 'CACCT', 'CACCC', 'CACCG', 'CACGA', 'CACGT', 'CACGC', 'CACGG', 'CAGAA', 'CAGAT', 'CAGAC', 'CAGAG', 'CAGTA', 'CAGTT', 'CAGTC', 'CAGTG', 'CAGCA', 'CAGCT', 'CAGCC', 'CAGCG', 'CAGGA', 'CAGGT', 'CAGGC', 'CAGGG', 'CTAAA', 'CTAAT', 'CTAAC', 'CTAAG', 'CTATA', 'CTATT', 'CTATC', 'CTATG', 'CTACA', 'CTACT', 'CTACC', 'CTACG', 'CTAGA', 'CTAGT', 'CTAGC', 'CTAGG', 'CTTAA', 'CTTAT', 'CTTAC', 'CTTAG', 'CTTTA', 'CTTTT', 'CTTTC', 'CTTTG', 'CTTCA', 'CTTCT', 'CTTCC', 'CTTCG', 'CTTGA', 'CTTGT', 'CTTGC', 'CTTGG', 'CTCAA', 'CTCAT', 'CTCAC', 'CTCAG', 'CTCTA', 'CTCTT', 'CTCTC', 'CTCTG', 'CTCCA', 'CTCCT', 'CTCCC', 'CTCCG', 'CTCGA', 'CTCGT', 'CTCGC', 'CTCGG', 'CTGAA', 'CTGAT', 'CTGAC', 'CTGAG', 'CTGTA', 'CTGTT', 'CTGTC', 'CTGTG', 'CTGCA', 'CTGCT', 'CTGCC', 'CTGCG', 'CTGGA', 'CTGGT', 'CTGGC', 'CTGGG', 'CCAAA', 'CCAAT', 'CCAAC', 'CCAAG', 'CCATA', 'CCATT', 'CCATC', 'CCATG', 'CCACA', 'CCACT', 'CCACC', 'CCACG', 'CCAGA', 'CCAGT', 'CCAGC', 'CCAGG', 'CCTAA', 'CCTAT', 'CCTAC', 'CCTAG', 'CCTTA', 'CCTTT', 'CCTTC', 'CCTTG', 'CCTCA', 'CCTCT', 'CCTCC', 'CCTCG', 'CCTGA', 'CCTGT', 'CCTGC', 'CCTGG', 'CCCAA', 'CCCAT', 'CCCAC', 'CCCAG', 'CCCTA', 'CCCTT', 'CCCTC', 'CCCTG', 'CCCCA', 'CCCCT', 'CCCCC', 'CCCCG', 'CCCGA', 'CCCGT', 'CCCGC', 'CCCGG', 'CCGAA', 'CCGAT', 'CCGAC', 'CCGAG', 'CCGTA', 'CCGTT', 'CCGTC', 'CCGTG', 'CCGCA', 'CCGCT', 'CCGCC', 'CCGCG', 'CCGGA', 'CCGGT', 'CCGGC', 'CCGGG', 'CGAAA', 'CGAAT', 'CGAAC', 'CGAAG', 'CGATA', 'CGATT', 'CGATC', 'CGATG', 'CGACA', 'CGACT', 'CGACC', 'CGACG', 'CGAGA', 'CGAGT', 'CGAGC', 'CGAGG', 'CGTAA', 'CGTAT', 'CGTAC', 'CGTAG', 'CGTTA', 'CGTTT', 'CGTTC', 'CGTTG', 'CGTCA', 'CGTCT', 'CGTCC', 'CGTCG', 'CGTGA', 'CGTGT', 'CGTGC', 'CGTGG', 'CGCAA', 'CGCAT', 'CGCAC', 'CGCAG', 'CGCTA', 'CGCTT', 'CGCTC', 'CGCTG', 'CGCCA', 'CGCCT', 'CGCCC', 'CGCCG', 'CGCGA', 'CGCGT', 'CGCGC', 'CGCGG', 'CGGAA', 'CGGAT', 'CGGAC', 'CGGAG', 'CGGTA', 'CGGTT', 'CGGTC', 'CGGTG', 'CGGCA', 'CGGCT', 'CGGCC', 'CGGCG', 'CGGGA', 'CGGGT', 'CGGGC', 'CGGGG', 'GAAAA', 'GAAAT', 'GAAAC', 'GAAAG', 'GAATA', 'GAATT', 'GAATC', 'GAATG', 'GAACA', 'GAACT', 'GAACC', 'GAACG', 'GAAGA', 'GAAGT', 'GAAGC', 'GAAGG', 'GATAA', 'GATAT', 'GATAC', 'GATAG', 'GATTA', 'GATTT', 'GATTC', 'GATTG', 'GATCA', 'GATCT', 'GATCC', 'GATCG', 'GATGA', 'GATGT', 'GATGC', 'GATGG', 'GACAA', 'GACAT', 'GACAC', 'GACAG', 'GACTA', 'GACTT', 'GACTC', 'GACTG', 'GACCA', 'GACCT', 'GACCC', 'GACCG', 'GACGA', 'GACGT', 'GACGC', 'GACGG', 'GAGAA', 'GAGAT', 'GAGAC', 'GAGAG', 'GAGTA', 'GAGTT', 'GAGTC', 'GAGTG', 'GAGCA', 'GAGCT', 'GAGCC', 'GAGCG', 'GAGGA', 'GAGGT', 'GAGGC', 'GAGGG', 'GTAAA', 'GTAAT', 'GTAAC', 'GTAAG', 'GTATA', 'GTATT', 'GTATC', 'GTATG', 'GTACA', 'GTACT', 'GTACC', 'GTACG', 'GTAGA', 'GTAGT', 'GTAGC', 'GTAGG', 'GTTAA', 'GTTAT', 'GTTAC', 'GTTAG', 'GTTTA', 'GTTTT', 'GTTTC', 'GTTTG', 'GTTCA', 'GTTCT', 'GTTCC', 'GTTCG', 'GTTGA', 'GTTGT', 'GTTGC', 'GTTGG', 'GTCAA', 'GTCAT', 'GTCAC', 'GTCAG', 'GTCTA', 'GTCTT', 'GTCTC', 'GTCTG', 'GTCCA', 'GTCCT', 'GTCCC', 'GTCCG', 'GTCGA', 'GTCGT', 'GTCGC', 'GTCGG', 'GTGAA', 'GTGAT', 'GTGAC', 'GTGAG', 'GTGTA', 'GTGTT', 'GTGTC', 'GTGTG', 'GTGCA', 'GTGCT', 'GTGCC', 'GTGCG', 'GTGGA', 'GTGGT', 'GTGGC', 'GTGGG', 'GCAAA', 'GCAAT', 'GCAAC', 'GCAAG', 'GCATA', 'GCATT', 'GCATC', 'GCATG', 'GCACA', 'GCACT', 'GCACC', 'GCACG', 'GCAGA', 'GCAGT', 'GCAGC', 'GCAGG', 'GCTAA', 'GCTAT', 'GCTAC', 'GCTAG', 'GCTTA', 'GCTTT', 'GCTTC', 'GCTTG', 'GCTCA', 'GCTCT', 'GCTCC', 'GCTCG', 'GCTGA', 'GCTGT', 'GCTGC', 'GCTGG', 'GCCAA', 'GCCAT', 'GCCAC', 'GCCAG', 'GCCTA', 'GCCTT', 'GCCTC', 'GCCTG', 'GCCCA', 'GCCCT', 'GCCCC', 'GCCCG', 'GCCGA', 'GCCGT', 'GCCGC', 'GCCGG', 'GCGAA', 'GCGAT', 'GCGAC', 'GCGAG', 'GCGTA', 'GCGTT', 'GCGTC', 'GCGTG', 'GCGCA', 'GCGCT', 'GCGCC', 'GCGCG', 'GCGGA', 'GCGGT', 'GCGGC', 'GCGGG', 'GGAAA', 'GGAAT', 'GGAAC', 'GGAAG', 'GGATA', 'GGATT', 'GGATC', 'GGATG', 'GGACA', 'GGACT', 'GGACC', 'GGACG', 'GGAGA', 'GGAGT', 'GGAGC', 'GGAGG', 'GGTAA', 'GGTAT', 'GGTAC', 'GGTAG', 'GGTTA', 'GGTTT', 'GGTTC', 'GGTTG', 'GGTCA', 'GGTCT', 'GGTCC', 'GGTCG', 'GGTGA', 'GGTGT', 'GGTGC', 'GGTGG', 'GGCAA', 'GGCAT', 'GGCAC', 'GGCAG', 'GGCTA', 'GGCTT', 'GGCTC', 'GGCTG', 'GGCCA', 'GGCCT', 'GGCCC', 'GGCCG', 'GGCGA', 'GGCGT', 'GGCGC', 'GGCGG', 'GGGAA', 'GGGAT', 'GGGAC', 'GGGAG', 'GGGTA', 'GGGTT', 'GGGTC', 'GGGTG', 'GGGCA', 'GGGCT', 'GGGCC', 'GGGCG', 'GGGGA', 'GGGGT', 'GGGGC', 'GGGGG')
i=-1
retro=0
DNA=0
gypsy=0
copia=0
LINE=0
SINE=0
L1=0
R1=0
with open("unknowns.fasta",'w') as pqr:
 with open("knowns.fasta",'w') as xyz:
  for name in fastanames1:
    k=fastanames1.index(name)
    entry1=fastalines1[k]
    entry2=fastalines1[k+1]
    fastaseq=""
    for m in range(entry1+1,entry2):
      fastaseq+=linecache.getline("itsp_" + "library_" + nameofg + ".fasta",m)
    if ("Unknown" not in name):
      xyz.write("%s%s\n%s\n" %(">",name,fastaseq))
    else:
      pqr.write("%s%s\n%s\n" %(">",name,fastaseq))
xyz.close()
pqr.close()

fastanames2=[]
fastalines2=[]
i=1
q=0
for row in open("knowns.fasta",'r'):
  if row.startswith(">") == True:
        fastanames2.append(row[1:])
        fastalines2.append(i)
        q+=1
  i+=1
fastalines2.append(i)
print "No. of known sequences:",q              

for name in fastanames2:
    k=fastanames2.index(name)
    entry1=fastalines2[k]
    entry2=fastalines2[k+1]
    i+=1
    dnaseq=""
    for m in range(entry1+1,entry2):
      dnaseq+=linecache.getline("knowns.fasta",m)
    #querynames.append(seq_record.id)
    if ("Gypsy" in name) or ("gypsy" in name) or ("GYPSY" in name) or ("PtAppalachian" in name) or ("PtTalladega" in name) or ("PtAngelina" in name) or ("PtOzark" in name) or ("PtOuachita" in name) or ("PtBastrop" in name):
       gypsy+=1
       retro+=1
    elif ("Copia" in name) or ("copia" in name) or ("COPIA" in name) or ("PtCumberland" in name) or ("PtConagree" in name) or ("PtPineywoods" in name):
       copia+=1
       retro+=1
    elif "LINE" in name:
       LINE+=1
       retro+=1
    elif ("retro" in name) or ("Retro" in name):
       retro+=1       
    elif "SINE" in name:
       SINE+=1
       retro+=1
    elif "LTR" in name:
       retro+=1
    elif "L1" in name:
       retro+=1
       LINE+=1
       L1+=1
    elif "R1" in name:
       retro+=1
       LINE+=1
       R1+=1
    elif "DNA" in name:
       DNA+=1
    elif ("Helitron" in name) or ("helitron" in name):
       DNA+=1
    elif ("Harbinger" in name) or ("harbinger" in name) or ("HARB" in name):
       DNA+=1
    elif ("EnSpm" in name):
       DNA+=1
    elif ("hAT" in name):
       DNA+=1
    elif ("MuDR" in name):
       DNA+=1
    elif ("Sola" in name) or ("sola" in name):
       DNA+=1   
    elif ("STOWAWAY" in name):
       DNA+=1
    elif ("Crypton" in name):
       DNA+=1
    else:
       DNA+=1
       print name
print "Total number of DNA elements:",DNA
print "Total number of retro elements:",retro
print "Gypsy:",gypsy," Copia:",copia," LINE:",LINE," SINE:",SINE," L1:",L1," R1:",R1
check=np.zeros(len(fastanames2))
dnano=0
ltrno=0
nonltrno=0

# Selects 2/3rd of known repeat sequences as query dataset
if DNA > retro:
   querycutoff=retro
if retro >= DNA:   
   querycutoff=DNA
with open("querydataset.fa",'w') as xyz:
    j=0
    m=-1
    for name in fastanames2:
       k=fastanames2.index(name)
       entry1=fastalines2[k]
       entry2=fastalines2[k+1]
       dnaseq=""
       m+=1
       for n in range(entry1+1,entry2):
           dnaseq+=linecache.getline("knowns.fasta",n)     
       if (("DNA" in name) or ("Helitron" in name) or ("helitron" in name) or ("Harbinger" in name) or ("harbinger" in name) or ("HARB" in name) or ("EnSpm" in name) or ("hAT" in name) or ("MuDR" in name) or ("Sola" in name) or ("sola" in name) or ("STOWAWAY" in name) or ("Crypton" in name)) and j<(querycutoff*2)/3: #
         check[m]=1
         j+=1
         dnano+=1
         xyz.write("%s%s\n%s\n" %(">",name,dnaseq))
    j=0
    s=0
    m=-1
    for name in fastanames2:
       k=fastanames2.index(name)
       entry1=fastalines2[k]
       entry2=fastalines2[k+1]
       dnaseq=""
       m+=1
       for n in range(entry1+1,entry2):
           dnaseq+=linecache.getline("knowns.fasta",n)     
       if (("Gypsy" in name) or ("gypsy" in name) or ("GYPSY" in name) or ("PtAppalachian" in name) or ("PtTalladega" in name) or ("PtAngelina" in name) or ("PtOzark" in name) or ("PtOuachita" in name) or ("PtBastrop" in name) or ("Copia" in name) or ("copia" in name) or ("COPIA" in name) or ("PtCumberland" in name) or ("PtConagree" in name) or ("PtPineywoods" in name) or ("LTR" in name)) and j<(querycutoff*2)/6:
         check[m]=1
         j+=1
         ltrno+=1
         xyz.write("%s%s\n%s\n" %(">",name,dnaseq))
       if (("LINE" in name) or ("SINE" in name) or ("L1" in name) or ("R1" in name)) and s<(querycutoff*2)/6:
         check[m]=1
         s+=1
         nonltrno+=1
         xyz.write("%s%s\n%s\n" %(">",name,dnaseq))
xyz.close()
print "Number of DNA elements in query dataset:",dnano
print "Number of LTR elements in query dataset:",ltrno
print "Number of non-LTR elements in query dataset:",s
countcheck=Counter(check)
print "Length of query dataset according to check:",countcheck[1]
with open("referencedataset.fa",'w') as xyz:
 j=0
 for m in check:
    if int(m)==0:
      name=fastanames2[j]
      entry1=fastalines2[j]
      entry2=fastalines2[j+1]
      dnaseq=""
      for n in range(entry1+1,entry2):
          dnaseq+=linecache.getline("knowns.fasta",n)      
      xyz.write("%s%s\n%s\n" %(">",name,dnaseq))
      j+=1
xyz.close()
print "No. of reference sequences:",j

fastanames3=[]
fastalines3=[]
i=1
for row in open("querydataset.fa",'r'): # Can give all_plantupdated.ref
  if row.startswith(">") == True:
        fastanames3.append(row[1:])
        fastalines3.append(i)
  i+=1
fastalines3.append(i)

fastanames4=[]
fastalines4=[]
i=1
for row in open("referencedataset.fa",'r'):
  if row.startswith(">") == True:
        fastanames4.append(row[1:])
        fastalines4.append(i)
  i+=1
fastalines4.append(i)

i=-1
print "Length of query dataset is:",len(fastanames3)
querydata=np.zeros((1024,len(fastanames3)))  # 4096 for hexanucleotides
querynames=[]  
for name in fastanames3:
    k=fastanames3.index(name)
    entry1=fastalines3[k]
    entry2=fastalines3[k+1]
    dnaseq=""
    for n in range(entry1+1,entry2):
       dnaseq+=linecache.getline("querydataset.fa",n)    
    i+=1
    if ("Gypsy" in name) or ("gypsy" in name) or ("GYPSY" in name):
       querynames.append("Retro")
    elif ("Copia" in name) or ("copia" in name) or ("COPIA" in name):
       querynames.append("Retro")
    elif "LINE" in name:
       querynames.append("Retro")
    elif ("PtCumberland" in name) or ("PtPineywoods" in name) or ("PtAppalachian" in name) or ("PtTalladega" in name) or ("PtAngelina" in name) or ("PtOzark" in name) or ("PtOuachita" in name) or ("PtBastrop" in name) or ("PtConagree" in name):
       querynames.append("Retro")       
    elif "SINE" in name:
       querynames.append("Retro")
    elif ("LTR" in name) or ("Retro" in name):
       querynames.append("Retro")
    elif "L1" in name:
       querynames.append("Retro")
    elif "R1" in name:
       querynames.append("Retro")
    elif "DNA" in name:
       querynames.append("DNA")
    elif ("Helitron" in name) or ("helitron" in name):
       querynames.append("DNA")
    elif ("Harbinger" in name) or ("harbinger" in name):
       querynames.append("DNA")
    elif ("EnSpm" in name):
       querynames.append("DNA")
    elif ("HARB" in name):
       querynames.append("DNA")       
    elif ("MuDR" in name):
       querynames.append("DNA")
    elif ("hAT" in name):
       querynames.append("DNA")
    elif ("Crypton" in name):
       querynames.append("DNA")       
    elif ("Sola" in name) or ("sola" in name):
       querynames.append("DNA")   
    elif ("STOWAWAY" in name):
       querynames.append("DNA")
    else:
       querynames.append("DNA")
       print "Exceptions:",name    
    j=-1
    for nuc in tetranucleotides:
        j+=1
        try:
          querydata[j,i]=dnaseq.count(nuc) # Calculating tetranucleotide frequency
        except IndexError:
          print j,i
querydata=np.transpose(querydata)
print "Query data is:",querydata
print "Query name is:",querynames

referencedata=np.zeros((1024,len(fastanames4))) # 4096 for hexanucleotides
referencenames=[]
i=-1
for name in fastanames4:
    k=fastanames4.index(name)
    entry1=fastalines4[k]
    entry2=fastalines4[k+1]
    dnaseq=""
    for n in range(entry1+1,entry2):
       dnaseq+=linecache.getline("referencedataset.fa",n) 
    i+=1
    if ("Gypsy" in name) or ("gypsy" in name) or ("GYPSY" in name):
       referencenames.append("Retro")
    elif ("Copia" in name) or ("copia" in name) or ("COPIA" in name):
       referencenames.append("Retro")
    elif "LINE" in name:
       referencenames.append("Retro")
    elif ("PtCumberland" in name) or ("PtPineywoods" in name) or ("PtAppalachian" in name) or ("PtTalladega" in name) or ("PtAngelina" in name) or ("PtOzark" in name) or ("PtOuachita" in name) or ("PtBastrop" in name) or ("PtConagree" in name):
       referencenames.append("Retro")       
    elif "SINE" in name:
       referencenames.append("Retro")
    elif ("LTR" in name) or ("Retro" in name):
       referencenames.append("Retro")
    elif "L1" in name:
       referencenames.append("Retro")
    elif "R1" in name:
       referencenames.append("Retro")
    elif "DNA" in name:
       referencenames.append("DNA")
    elif ("Helitron" in name) or ("helitron" in name):
       referencenames.append("DNA")
    elif ("Harbinger" in name) or ("harbinger" in name):
       referencenames.append("DNA")
    elif ("EnSpm" in name):
       referencenames.append("DNA")
    elif ("HARB" in name):
       referencenames.append("DNA")       
    elif ("MuDR" in name):
       referencenames.append("DNA")
    elif ("hAT" in name):
       referencenames.append("DNA")
    elif ("Crypton" in name):
       referencenames.append("DNA")       
    elif ("Sola" in name) or ("sola" in name):
       referencenames.append("DNA")
    else:
       referencenames.append("DNA")
       print "Exceptions:",name    
    j=-1
    for nuc in tetranucleotides:
        j+=1
        referencedata[j,i]=dnaseq.count(nuc) # Calculating tetranucleotide frequency
referencedata=np.transpose(referencedata)
print "Reference data is:",referencedata

# Validating datasets

clf = svm.SVC(gamma=0.0070, C=5)
clf.fit(querydata, querynames)
referencepredict=clf.predict(referencedata)
print "Confusion matrix for DNA/LTR retrotransposon is:",confusion_matrix(referencenames, referencepredict)
print "Accuracy for DNA/LTR retrotransposon is:", accuracy_score(referencenames, referencepredict)

# Starting predictions of unknowns

fastanames4=[]
fastalines4=[]
i=1
for row in open("unknowns.fasta",'r'):
  if row.startswith(">") == True:
        fastanames4.append(row[1:])
        fastalines4.append(i)
  i+=1
fastalines4.append(i)

referencedata=np.zeros((1024,len(fastanames4))) # 4096 for hexanucleotides
i=-1
for name in fastanames4:
    k=fastanames4.index(name)
    entry1=fastalines4[k]
    entry2=fastalines4[k+1]
    dnaseq=""
    for n in range(entry1+1,entry2):
       dnaseq+=linecache.getline("unknowns.fasta",n) 
    i+=1
    for nuc in tetranucleotides:
        j+=1
        referencedata[j,i]=dnaseq.count(nuc) # Calculating tetranucleotide frequency
referencedata=np.transpose(referencedata)
print "Reference data is:",referencedata
# Machine learning starts

# Predicting DNA transposons and retrotransposons 

clf = svm.SVC(gamma=0.0070, C=5)
clf.fit(querydata, querynames)
referencepredict=clf.predict(referencedata)
repeatclassification=[[] for i in range(len(fastanames4))]
retros=[]
i=0
repmods=0  #Use this variable for adding ML classifications to 'repeatclassification' list
for repclass in referencepredict:
    if repclass=="DNA":
       repeatclassification[repmods]=[fastanames4[i],"DNA"]
       repmods+=1
    if repclass=="Retro":
       retros.append(fastanames4[i])     
    i+=1   
# Removing temporary files
#os.remove("itsp_" + "library_" + nameofg + ".fasta")
