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

time_start = time.time()
parser = argparse.ArgumentParser(
     prog='repeatpipeline',
     usage='''python repeatpipeline.py --fasta [Name of fasta file] --name [name of genome (optional)] --path [Path to fasta file (optional)] --cpus [Number of CPU's (optional)] --repbase [Name of reference database (optional)] --pog [Percentage of genome for de novo repeat detection (optional)]''',
     description='''This program is a complete pipeline for identifying and characterizing transposable elements in a given genome.''',
     epilog='''It requires numpy, matplotlib, scikit and biopython libraries for execution.''')
parser.add_argument('--fasta', type=str, help='The query fasta file', required=True)
parser.add_argument('--name', type=str, help='The name of genome (optional)', required=False)
parser.add_argument('--pog', type=int, help='Percentage of genome for de novo repeat detection, please give whole numbers (optional)', required=False)
parser.add_argument('--path', type=str, help='Path to input files [Not needed if genome file in same directory as this script]', required=False)
parser.add_argument('--cpus', type=int, help="Number of CPU's (optional)", required=False)
parser.add_argument('--repbase', type=str, help='Reference database (from repbase) (optional)', required=False)
args = parser.parse_args()

filename1=args.fasta
path=args.path
nameofg=args.name
cpus=args.cpus
pog=args.pog
repbasefile=args.repbase
if path == None:
   filename5=open(filename1,"r")
else:
   filename5=open(os.path.join(path, filename1), "r")

if nameofg == None:
   nameofg="temp"

if cpus == None:
   cpus=1    

# Making a new directory
status, homedir = commands.getstatusoutput("pwd")
foldername="repeat_" + nameofg + "_" + str(datetime.date.today()) + "_"  + str(time.time())
print "Name of folder:",foldername
os.system("mkdir " + foldername)
os.chdir(foldername)

# Running repeatmodeler
os.system("mkdir repeatmodeler")
os.chdir("repeatmodeler")
if pog == None:
   if path == None:
        os.system("BuildDatabase -name " + nameofg + " -engine ncbi " + filename1)
   else:
        os.system("BuildDatabase -name " + nameofg + " -engine ncbi " + path + "/" + filename1)
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

   fastanames=[]
   fastalines=[]
   i=1
   for row in filename8:
     if row.startswith(">") == True:
        fastanames.append(row[1:])
        fastalines.append(i)
     i+=1
   fastalines.append(i)
   # Calulate total number of bases

   genome=""

   for name in fastanames:
       k=fastanames.index(name)
       entry1=fastalines[k]
       entry2=fastalines[k+1]
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
         while (breaks < splice) and j < len(fastanames):
             entry1=fastalines[j]
             entry2=fastalines[j+1]
             print entry1,entry2
             sequence=""
             for m in range(entry1+1,entry2):
                 sequence+=linecache.getline(filename9,m)
             ijk.write("%s%s\n%s\n" % (">",fastanames[j].rstrip('\r\n'), sequence.rstrip('\r\n'))) # Making the fasta file
             #print j,fastanames[j]
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
   print "BuildDatabase -name " + nameofg + " -engine ncbi repeatmodelertrainingset.fasta"
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
    os.system("cp ../../../all_plantupdated.ref all_plantupdated.ref")
    repbasefile="all_plantupdated.ref"
else:
    os.system("cp ../../../" + repbasefile + " " + repbasefile)
os.system("censor.ncbi consensi.fa.classified -lib " + repbasefile  + " -tab -no_simple -s")

# Incorporating CENSOR results into the main repeatmodeler library

with open("unknown8080.txt",'w') as xyz:
 for line in open("consensi.fa.classified.map",'r'):
    if "Unknown" in line:
       line2=line.split()
       if float(line2[7]) > 0.8 and (int(line2[2])-int(line2[1]) >= 80): # Applying the 80/80 cut-off.
          xyz.write("%s" %(line))
xyz.close()

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

# Running repeatmasker
os.system("mkdir repeatmasker")
os.chdir("repeatmasker")
status, repeatmaskerdir = commands.getstatusoutput("pwd")
repeatmaskerstart=time.time()
if path == None:
   os.system("RepeatMasker " + str(homedir) + "/" + str(filename1) + " -dir " + str(repeatmaskerdir) + " -lib " + str(repeatmodelerdir) + "/" + str(nameofg) + "_completelib.fasta" + " -pa " + str(cpus) + " -gff -a -noisy -low")
else:   
   os.system("RepeatMasker " + str(path) + "/" + str(filename1) + " -dir " + str(repeatmaskerdir) + " -lib " + str(repeatmodelerdir) + "/" + str(nameofg) + "_completelib.fasta" + " -pa " + str(cpus) + " -gff -a -noisy -low")
print "Time taken to run RepeatMasker is: ",time.time() - repeatmaskerstart
