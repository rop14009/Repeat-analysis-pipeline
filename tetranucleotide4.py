# python tetranucleotide4.py --query Pila2.5f_consensi.fa.classified
from Bio import SeqIO
import argparse
import numpy as np
import linecache
import time
from sklearn import svm
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score

time_start = time.time()
parser = argparse.ArgumentParser(
     prog='tetranucleotide',
     usage='''python tetranucleotide4.py --query [queryfile] --output [outputfile]''',
     description='''This program calculates tetranucleotide frequencies of the given DNA sequences.''',
     epilog='''It requires numpy, matplotlib, scikit and biopython libraries for execution''')
parser.add_argument('--query', type=str, help='The query fasta file', required=True)
parser.add_argument('--output', type=str, help='The output file (optional)', required=False)
args = parser.parse_args()

filename1=args.query
filename2=args.output

seq_records1 = SeqIO.parse(filename1, "fasta")
seq_records1 = list(seq_records1)

tetranucleotides =('AAAAA', 'AAAAT', 'AAAAC', 'AAAAG', 'AAATA', 'AAATT', 'AAATC', 'AAATG', 'AAACA', 'AAACT', 'AAACC', 'AAACG', 'AAAGA', 'AAAGT', 'AAAGC', 'AAAGG', 'AATAA', 'AATAT', 'AATAC', 'AATAG', 'AATTA', 'AATTT', 'AATTC', 'AATTG', 'AATCA', 'AATCT', 'AATCC', 'AATCG', 'AATGA', 'AATGT', 'AATGC', 'AATGG', 'AACAA', 'AACAT', 'AACAC', 'AACAG', 'AACTA', 'AACTT', 'AACTC', 'AACTG', 'AACCA', 'AACCT', 'AACCC', 'AACCG', 'AACGA', 'AACGT', 'AACGC', 'AACGG', 'AAGAA', 'AAGAT', 'AAGAC', 'AAGAG', 'AAGTA', 'AAGTT', 'AAGTC', 'AAGTG', 'AAGCA', 'AAGCT', 'AAGCC', 'AAGCG', 'AAGGA', 'AAGGT', 'AAGGC', 'AAGGG', 'ATAAA', 'ATAAT', 'ATAAC', 'ATAAG', 'ATATA', 'ATATT', 'ATATC', 'ATATG', 'ATACA', 'ATACT', 'ATACC', 'ATACG', 'ATAGA', 'ATAGT', 'ATAGC', 'ATAGG', 'ATTAA', 'ATTAT', 'ATTAC', 'ATTAG', 'ATTTA', 'ATTTT', 'ATTTC', 'ATTTG', 'ATTCA', 'ATTCT', 'ATTCC', 'ATTCG', 'ATTGA', 'ATTGT', 'ATTGC', 'ATTGG', 'ATCAA', 'ATCAT', 'ATCAC', 'ATCAG', 'ATCTA', 'ATCTT', 'ATCTC', 'ATCTG', 'ATCCA', 'ATCCT', 'ATCCC', 'ATCCG', 'ATCGA', 'ATCGT', 'ATCGC', 'ATCGG', 'ATGAA', 'ATGAT', 'ATGAC', 'ATGAG', 'ATGTA', 'ATGTT', 'ATGTC', 'ATGTG', 'ATGCA', 'ATGCT', 'ATGCC', 'ATGCG', 'ATGGA', 'ATGGT', 'ATGGC', 'ATGGG', 'ACAAA', 'ACAAT', 'ACAAC', 'ACAAG', 'ACATA', 'ACATT', 'ACATC', 'ACATG', 'ACACA', 'ACACT', 'ACACC', 'ACACG', 'ACAGA', 'ACAGT', 'ACAGC', 'ACAGG', 'ACTAA', 'ACTAT', 'ACTAC', 'ACTAG', 'ACTTA', 'ACTTT', 'ACTTC', 'ACTTG', 'ACTCA', 'ACTCT', 'ACTCC', 'ACTCG', 'ACTGA', 'ACTGT', 'ACTGC', 'ACTGG', 'ACCAA', 'ACCAT', 'ACCAC', 'ACCAG', 'ACCTA', 'ACCTT', 'ACCTC', 'ACCTG', 'ACCCA', 'ACCCT', 'ACCCC', 'ACCCG', 'ACCGA', 'ACCGT', 'ACCGC', 'ACCGG', 'ACGAA', 'ACGAT', 'ACGAC', 'ACGAG', 'ACGTA', 'ACGTT', 'ACGTC', 'ACGTG', 'ACGCA', 'ACGCT', 'ACGCC', 'ACGCG', 'ACGGA', 'ACGGT', 'ACGGC', 'ACGGG', 'AGAAA', 'AGAAT', 'AGAAC', 'AGAAG', 'AGATA', 'AGATT', 'AGATC', 'AGATG', 'AGACA', 'AGACT', 'AGACC', 'AGACG', 'AGAGA', 'AGAGT', 'AGAGC', 'AGAGG', 'AGTAA', 'AGTAT', 'AGTAC', 'AGTAG', 'AGTTA', 'AGTTT', 'AGTTC', 'AGTTG', 'AGTCA', 'AGTCT', 'AGTCC', 'AGTCG', 'AGTGA', 'AGTGT', 'AGTGC', 'AGTGG', 'AGCAA', 'AGCAT', 'AGCAC', 'AGCAG', 'AGCTA', 'AGCTT', 'AGCTC', 'AGCTG', 'AGCCA', 'AGCCT', 'AGCCC', 'AGCCG', 'AGCGA', 'AGCGT', 'AGCGC', 'AGCGG', 'AGGAA', 'AGGAT', 'AGGAC', 'AGGAG', 'AGGTA', 'AGGTT', 'AGGTC', 'AGGTG', 'AGGCA', 'AGGCT', 'AGGCC', 'AGGCG', 'AGGGA', 'AGGGT', 'AGGGC', 'AGGGG', 'TAAAA', 'TAAAT', 'TAAAC', 'TAAAG', 'TAATA', 'TAATT', 'TAATC', 'TAATG', 'TAACA', 'TAACT', 'TAACC', 'TAACG', 'TAAGA', 'TAAGT', 'TAAGC', 'TAAGG', 'TATAA', 'TATAT', 'TATAC', 'TATAG', 'TATTA', 'TATTT', 'TATTC', 'TATTG', 'TATCA', 'TATCT', 'TATCC', 'TATCG', 'TATGA', 'TATGT', 'TATGC', 'TATGG', 'TACAA', 'TACAT', 'TACAC', 'TACAG', 'TACTA', 'TACTT', 'TACTC', 'TACTG', 'TACCA', 'TACCT', 'TACCC', 'TACCG', 'TACGA', 'TACGT', 'TACGC', 'TACGG', 'TAGAA', 'TAGAT', 'TAGAC', 'TAGAG', 'TAGTA', 'TAGTT', 'TAGTC', 'TAGTG', 'TAGCA', 'TAGCT', 'TAGCC', 'TAGCG', 'TAGGA', 'TAGGT', 'TAGGC', 'TAGGG', 'TTAAA', 'TTAAT', 'TTAAC', 'TTAAG', 'TTATA', 'TTATT', 'TTATC', 'TTATG', 'TTACA', 'TTACT', 'TTACC', 'TTACG', 'TTAGA', 'TTAGT', 'TTAGC', 'TTAGG', 'TTTAA', 'TTTAT', 'TTTAC', 'TTTAG', 'TTTTA', 'TTTTT', 'TTTTC', 'TTTTG', 'TTTCA', 'TTTCT', 'TTTCC', 'TTTCG', 'TTTGA', 'TTTGT', 'TTTGC', 'TTTGG', 'TTCAA', 'TTCAT', 'TTCAC', 'TTCAG', 'TTCTA', 'TTCTT', 'TTCTC', 'TTCTG', 'TTCCA', 'TTCCT', 'TTCCC', 'TTCCG', 'TTCGA', 'TTCGT', 'TTCGC', 'TTCGG', 'TTGAA', 'TTGAT', 'TTGAC', 'TTGAG', 'TTGTA', 'TTGTT', 'TTGTC', 'TTGTG', 'TTGCA', 'TTGCT', 'TTGCC', 'TTGCG', 'TTGGA', 'TTGGT', 'TTGGC', 'TTGGG', 'TCAAA', 'TCAAT', 'TCAAC', 'TCAAG', 'TCATA', 'TCATT', 'TCATC', 'TCATG', 'TCACA', 'TCACT', 'TCACC', 'TCACG', 'TCAGA', 'TCAGT', 'TCAGC', 'TCAGG', 'TCTAA', 'TCTAT', 'TCTAC', 'TCTAG', 'TCTTA', 'TCTTT', 'TCTTC', 'TCTTG', 'TCTCA', 'TCTCT', 'TCTCC', 'TCTCG', 'TCTGA', 'TCTGT', 'TCTGC', 'TCTGG', 'TCCAA', 'TCCAT', 'TCCAC', 'TCCAG', 'TCCTA', 'TCCTT', 'TCCTC', 'TCCTG', 'TCCCA', 'TCCCT', 'TCCCC', 'TCCCG', 'TCCGA', 'TCCGT', 'TCCGC', 'TCCGG', 'TCGAA', 'TCGAT', 'TCGAC', 'TCGAG', 'TCGTA', 'TCGTT', 'TCGTC', 'TCGTG', 'TCGCA', 'TCGCT', 'TCGCC', 'TCGCG', 'TCGGA', 'TCGGT', 'TCGGC', 'TCGGG', 'TGAAA', 'TGAAT', 'TGAAC', 'TGAAG', 'TGATA', 'TGATT', 'TGATC', 'TGATG', 'TGACA', 'TGACT', 'TGACC', 'TGACG', 'TGAGA', 'TGAGT', 'TGAGC', 'TGAGG', 'TGTAA', 'TGTAT', 'TGTAC', 'TGTAG', 'TGTTA', 'TGTTT', 'TGTTC', 'TGTTG', 'TGTCA', 'TGTCT', 'TGTCC', 'TGTCG', 'TGTGA', 'TGTGT', 'TGTGC', 'TGTGG', 'TGCAA', 'TGCAT', 'TGCAC', 'TGCAG', 'TGCTA', 'TGCTT', 'TGCTC', 'TGCTG', 'TGCCA', 'TGCCT', 'TGCCC', 'TGCCG', 'TGCGA', 'TGCGT', 'TGCGC', 'TGCGG', 'TGGAA', 'TGGAT', 'TGGAC', 'TGGAG', 'TGGTA', 'TGGTT', 'TGGTC', 'TGGTG', 'TGGCA', 'TGGCT', 'TGGCC', 'TGGCG', 'TGGGA', 'TGGGT', 'TGGGC', 'TGGGG', 'CAAAA', 'CAAAT', 'CAAAC', 'CAAAG', 'CAATA', 'CAATT', 'CAATC', 'CAATG', 'CAACA', 'CAACT', 'CAACC', 'CAACG', 'CAAGA', 'CAAGT', 'CAAGC', 'CAAGG', 'CATAA', 'CATAT', 'CATAC', 'CATAG', 'CATTA', 'CATTT', 'CATTC', 'CATTG', 'CATCA', 'CATCT', 'CATCC', 'CATCG', 'CATGA', 'CATGT', 'CATGC', 'CATGG', 'CACAA', 'CACAT', 'CACAC', 'CACAG', 'CACTA', 'CACTT', 'CACTC', 'CACTG', 'CACCA', 'CACCT', 'CACCC', 'CACCG', 'CACGA', 'CACGT', 'CACGC', 'CACGG', 'CAGAA', 'CAGAT', 'CAGAC', 'CAGAG', 'CAGTA', 'CAGTT', 'CAGTC', 'CAGTG', 'CAGCA', 'CAGCT', 'CAGCC', 'CAGCG', 'CAGGA', 'CAGGT', 'CAGGC', 'CAGGG', 'CTAAA', 'CTAAT', 'CTAAC', 'CTAAG', 'CTATA', 'CTATT', 'CTATC', 'CTATG', 'CTACA', 'CTACT', 'CTACC', 'CTACG', 'CTAGA', 'CTAGT', 'CTAGC', 'CTAGG', 'CTTAA', 'CTTAT', 'CTTAC', 'CTTAG', 'CTTTA', 'CTTTT', 'CTTTC', 'CTTTG', 'CTTCA', 'CTTCT', 'CTTCC', 'CTTCG', 'CTTGA', 'CTTGT', 'CTTGC', 'CTTGG', 'CTCAA', 'CTCAT', 'CTCAC', 'CTCAG', 'CTCTA', 'CTCTT', 'CTCTC', 'CTCTG', 'CTCCA', 'CTCCT', 'CTCCC', 'CTCCG', 'CTCGA', 'CTCGT', 'CTCGC', 'CTCGG', 'CTGAA', 'CTGAT', 'CTGAC', 'CTGAG', 'CTGTA', 'CTGTT', 'CTGTC', 'CTGTG', 'CTGCA', 'CTGCT', 'CTGCC', 'CTGCG', 'CTGGA', 'CTGGT', 'CTGGC', 'CTGGG', 'CCAAA', 'CCAAT', 'CCAAC', 'CCAAG', 'CCATA', 'CCATT', 'CCATC', 'CCATG', 'CCACA', 'CCACT', 'CCACC', 'CCACG', 'CCAGA', 'CCAGT', 'CCAGC', 'CCAGG', 'CCTAA', 'CCTAT', 'CCTAC', 'CCTAG', 'CCTTA', 'CCTTT', 'CCTTC', 'CCTTG', 'CCTCA', 'CCTCT', 'CCTCC', 'CCTCG', 'CCTGA', 'CCTGT', 'CCTGC', 'CCTGG', 'CCCAA', 'CCCAT', 'CCCAC', 'CCCAG', 'CCCTA', 'CCCTT', 'CCCTC', 'CCCTG', 'CCCCA', 'CCCCT', 'CCCCC', 'CCCCG', 'CCCGA', 'CCCGT', 'CCCGC', 'CCCGG', 'CCGAA', 'CCGAT', 'CCGAC', 'CCGAG', 'CCGTA', 'CCGTT', 'CCGTC', 'CCGTG', 'CCGCA', 'CCGCT', 'CCGCC', 'CCGCG', 'CCGGA', 'CCGGT', 'CCGGC', 'CCGGG', 'CGAAA', 'CGAAT', 'CGAAC', 'CGAAG', 'CGATA', 'CGATT', 'CGATC', 'CGATG', 'CGACA', 'CGACT', 'CGACC', 'CGACG', 'CGAGA', 'CGAGT', 'CGAGC', 'CGAGG', 'CGTAA', 'CGTAT', 'CGTAC', 'CGTAG', 'CGTTA', 'CGTTT', 'CGTTC', 'CGTTG', 'CGTCA', 'CGTCT', 'CGTCC', 'CGTCG', 'CGTGA', 'CGTGT', 'CGTGC', 'CGTGG', 'CGCAA', 'CGCAT', 'CGCAC', 'CGCAG', 'CGCTA', 'CGCTT', 'CGCTC', 'CGCTG', 'CGCCA', 'CGCCT', 'CGCCC', 'CGCCG', 'CGCGA', 'CGCGT', 'CGCGC', 'CGCGG', 'CGGAA', 'CGGAT', 'CGGAC', 'CGGAG', 'CGGTA', 'CGGTT', 'CGGTC', 'CGGTG', 'CGGCA', 'CGGCT', 'CGGCC', 'CGGCG', 'CGGGA', 'CGGGT', 'CGGGC', 'CGGGG', 'GAAAA', 'GAAAT', 'GAAAC', 'GAAAG', 'GAATA', 'GAATT', 'GAATC', 'GAATG', 'GAACA', 'GAACT', 'GAACC', 'GAACG', 'GAAGA', 'GAAGT', 'GAAGC', 'GAAGG', 'GATAA', 'GATAT', 'GATAC', 'GATAG', 'GATTA', 'GATTT', 'GATTC', 'GATTG', 'GATCA', 'GATCT', 'GATCC', 'GATCG', 'GATGA', 'GATGT', 'GATGC', 'GATGG', 'GACAA', 'GACAT', 'GACAC', 'GACAG', 'GACTA', 'GACTT', 'GACTC', 'GACTG', 'GACCA', 'GACCT', 'GACCC', 'GACCG', 'GACGA', 'GACGT', 'GACGC', 'GACGG', 'GAGAA', 'GAGAT', 'GAGAC', 'GAGAG', 'GAGTA', 'GAGTT', 'GAGTC', 'GAGTG', 'GAGCA', 'GAGCT', 'GAGCC', 'GAGCG', 'GAGGA', 'GAGGT', 'GAGGC', 'GAGGG', 'GTAAA', 'GTAAT', 'GTAAC', 'GTAAG', 'GTATA', 'GTATT', 'GTATC', 'GTATG', 'GTACA', 'GTACT', 'GTACC', 'GTACG', 'GTAGA', 'GTAGT', 'GTAGC', 'GTAGG', 'GTTAA', 'GTTAT', 'GTTAC', 'GTTAG', 'GTTTA', 'GTTTT', 'GTTTC', 'GTTTG', 'GTTCA', 'GTTCT', 'GTTCC', 'GTTCG', 'GTTGA', 'GTTGT', 'GTTGC', 'GTTGG', 'GTCAA', 'GTCAT', 'GTCAC', 'GTCAG', 'GTCTA', 'GTCTT', 'GTCTC', 'GTCTG', 'GTCCA', 'GTCCT', 'GTCCC', 'GTCCG', 'GTCGA', 'GTCGT', 'GTCGC', 'GTCGG', 'GTGAA', 'GTGAT', 'GTGAC', 'GTGAG', 'GTGTA', 'GTGTT', 'GTGTC', 'GTGTG', 'GTGCA', 'GTGCT', 'GTGCC', 'GTGCG', 'GTGGA', 'GTGGT', 'GTGGC', 'GTGGG', 'GCAAA', 'GCAAT', 'GCAAC', 'GCAAG', 'GCATA', 'GCATT', 'GCATC', 'GCATG', 'GCACA', 'GCACT', 'GCACC', 'GCACG', 'GCAGA', 'GCAGT', 'GCAGC', 'GCAGG', 'GCTAA', 'GCTAT', 'GCTAC', 'GCTAG', 'GCTTA', 'GCTTT', 'GCTTC', 'GCTTG', 'GCTCA', 'GCTCT', 'GCTCC', 'GCTCG', 'GCTGA', 'GCTGT', 'GCTGC', 'GCTGG', 'GCCAA', 'GCCAT', 'GCCAC', 'GCCAG', 'GCCTA', 'GCCTT', 'GCCTC', 'GCCTG', 'GCCCA', 'GCCCT', 'GCCCC', 'GCCCG', 'GCCGA', 'GCCGT', 'GCCGC', 'GCCGG', 'GCGAA', 'GCGAT', 'GCGAC', 'GCGAG', 'GCGTA', 'GCGTT', 'GCGTC', 'GCGTG', 'GCGCA', 'GCGCT', 'GCGCC', 'GCGCG', 'GCGGA', 'GCGGT', 'GCGGC', 'GCGGG', 'GGAAA', 'GGAAT', 'GGAAC', 'GGAAG', 'GGATA', 'GGATT', 'GGATC', 'GGATG', 'GGACA', 'GGACT', 'GGACC', 'GGACG', 'GGAGA', 'GGAGT', 'GGAGC', 'GGAGG', 'GGTAA', 'GGTAT', 'GGTAC', 'GGTAG', 'GGTTA', 'GGTTT', 'GGTTC', 'GGTTG', 'GGTCA', 'GGTCT', 'GGTCC', 'GGTCG', 'GGTGA', 'GGTGT', 'GGTGC', 'GGTGG', 'GGCAA', 'GGCAT', 'GGCAC', 'GGCAG', 'GGCTA', 'GGCTT', 'GGCTC', 'GGCTG', 'GGCCA', 'GGCCT', 'GGCCC', 'GGCCG', 'GGCGA', 'GGCGT', 'GGCGC', 'GGCGG', 'GGGAA', 'GGGAT', 'GGGAC', 'GGGAG', 'GGGTA', 'GGGTT', 'GGGTC', 'GGGTG', 'GGGCA', 'GGGCT', 'GGGCC', 'GGGCG', 'GGGGA', 'GGGGT', 'GGGGC', 'GGGGG')
i=-1
retro=0
DNA=0
simple=0
RNA=0
gypsy=0
copia=0
LINE=0
SINE=0
L1=0
R1=0
with open("unknowns.fasta",'w') as pqr:
 with open("knowns.fasta",'w') as xyz:
  for seq_record in seq_records1:
    if ("Unknown" not in seq_record.id):
      xyz.write("%s%s\n%s\n" %(">",str(seq_record.id),str(seq_record.seq)))
    else:
      pqr.write("%s%s\n%s\n" %(">",str(seq_record.id),str(seq_record.seq)))
xyz.close()
pqr.close()
seq_records2 = SeqIO.parse("knowns.fasta", "fasta")
seq_records2 = list(seq_records2)

for seq_record in seq_records2:
    i+=1
    dnaseq = str(seq_record.seq)
    #querynames.append(seq_record.id)
    if ("Gypsy" in seq_record.id) or ("gypsy" in seq_record.id):
       gypsy+=1
       retro+=1
    elif ("Copia" in seq_record.id) or ("copia" in seq_record.id):
       copia+=1
       retro+=1
    elif "LINE" in seq_record.id:
       LINE+=1
       retro+=1
    elif "SINE" in seq_record.id:
       SINE+=1
       retro+=1
    elif "LTR" in seq_record.id:
       retro+=1
    elif "L1" in seq_record.id:
       retro+=1
       L1+=1
    elif "R1" in seq_record.id:
       retro+=1
       R1+=1
    elif "DNA" in seq_record.id:
       DNA+=1
    elif ("Helitron" in seq_record.id) or ("helitron" in seq_record.id):
       DNA+=1
    elif ("Harbinger" in seq_record.id) or ("harbinger" in seq_record.id):
       DNA+=1
    elif ("EnSpm" in seq_record.id):
       DNA+=1
    elif ("MuDR" in seq_record.id):
       DNA+=1
    elif ("Sola" in seq_record.id) or ("sola" in seq_record.id):
       DNA+=1   
    elif ("STOWAWAY" in seq_record.id):
       DNA+=1
    elif ("Satellite" in seq_record.id):
       simple+=1
    elif ("Simple_repeat" in seq_record.id):
       simple+=1
    elif ("rRNA" in seq_record.id):
       RNA+=1
    else:
       DNA+=1
       #print seq_record.id
print "No. of DNA elements:",DNA


check=np.zeros(len(seq_records2))
with open("querydataset.fa",'w') as xyz:
 if DNA < retro:
    j=0
    m=0
    for seq_record in seq_records2:
       if "DNA" in str(seq_record.id) and j<DNA/2.5: #
         check[m]=1
         j+=1
         m+=1
         xyz.write("%s%s\n%s\n" %(">",str(seq_record.id),str(seq_record.seq)))
    m=0
    for seq_record in seq_records2:
       if "Simple_repeat" in str(seq_record.id):
         check[m]=1
         m+=1
         j+=1
         continue
         #xyz.write("%s%s\n%s\n" %(">",str(seq_record.id),str(seq_record.seq)))
    m=0
    for seq_record in seq_records2:
       if "Satellite" in str(seq_record.id):
         check[m]=1
         m+=1
         j+=1
         continue
         #xyz.write("%s%s\n%s\n" %(">",str(seq_record.id),str(seq_record.seq)))
    m=0
    for seq_record in seq_records2:
       if "rRNA" in str(seq_record.id):
         check[m]=1
         m+=1
         j+=1
         continue
         #xyz.write("%s%s\n%s\n" %(">",str(seq_record.id),str(seq_record.seq)))
    m=0
    for seq_record in seq_records2:
       if (("Copia" in str(seq_record.id)) or ("Gypsy" in str(seq_record.id))) and j<DNA:
         check[m]=1
         m+=1
         j+=1
         xyz.write("%s%s\n%s\n" %(">",str(seq_record.id),str(seq_record.seq)))
xyz.close()
with open("referencedataset.fa",'w') as xyz:
 j=0
 for m in check:
    if int(m)==0:
     if ("Simple_repeat" in str(seq_records2[j].id)) or ("Satellite" in str(seq_records2[j].id)) or ("rRNA" in str(seq_records2[j].id)):
      continue
     else:
      xyz.write("%s%s\n%s\n" %(">",str(seq_records2[j].id),str(seq_records2[j].seq)))
      j+=1
xyz.close()
seq_records3 = SeqIO.parse("querydataset.fa", "fasta")
#seq_records3 = SeqIO.parse("all_plant.ref", "fasta")
seq_records3 = list(seq_records3)

seq_records4 = SeqIO.parse("referencedataset.fa", "fasta")
#seq_records4 = SeqIO.parse("knowns.fasta", "fasta")
seq_records4 = list(seq_records4)
i=-1
print "Length of query dataset is:",len(seq_records3)
querydata=np.zeros((1024,len(seq_records3)))  # 4096 for hexanucleotides
querynames=[]  
for seq_record in seq_records3:
    i+=1
    dnaseq = str(seq_record.seq)
    #querynames.append(seq_record.id)
    if ("Gypsy" in seq_record.id) or ("gypsy" in seq_record.id):
       querynames.append("Retro")
    elif ("Copia" in seq_record.id) or ("copia" in seq_record.id):
       querynames.append("Retro")
    elif "LINE" in seq_record.id:
       querynames.append("Retro")
    elif "SINE" in seq_record.id:
       querynames.append("Retro")
    elif ("LTR" in seq_record.id) or ("Retro" in seq_record.id):
       querynames.append("Retro")
    elif "L1" in seq_record.id:
       querynames.append("Retro")
    elif "R1" in seq_record.id:
       querynames.append("Retro")
    elif "DNA" in seq_record.id:
       querynames.append("DNA")
    elif ("Helitron" in seq_record.id) or ("helitron" in seq_record.id):
       querynames.append("DNA")
    elif ("Harbinger" in seq_record.id) or ("harbinger" in seq_record.id):
       querynames.append("DNA")
    elif ("EnSpm" in seq_record.id):
       querynames.append("DNA")
    elif ("MuDR" in seq_record.id):
       querynames.append("DNA")
    elif ("Sola" in seq_record.id) or ("sola" in seq_record.id):
       querynames.append("DNA")   
    elif ("STOWAWAY" in seq_record.id):
       querynames.append("DNA")
    elif ("Satellite" in seq_record.id):
       continue
    elif ("Simple_repeat" in seq_record.id):
       continue
    elif ("rRNA" in seq_record.id):
       continue
    else:
       querynames.append("DNA")
       print seq_record.id    
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

referencedata=np.zeros((1024,len(seq_records4))) # 4096 for hexanucleotides
referencenames=[]
i=-1
for seq_record in seq_records4:
    i+=1
    dnaseq = str(seq_record.seq)
    #referencenames.append(seq_record.id)
    if ("Gypsy" in seq_record.id) or ("gypsy" in seq_record.id):
       referencenames.append("Retro")
    elif ("Copia" in seq_record.id) or ("copia" in seq_record.id):
       referencenames.append("Retro")
    elif "LINE" in seq_record.id:
       referencenames.append("Retro")
    elif "SINE" in seq_record.id:
       referencenames.append("Retro")
    elif ("LTR" in seq_record.id) or ("Retro" in seq_record.id):
       referencenames.append("Retro")
    elif "L1" in seq_record.id:
       referencenames.append("Retro")
    elif "R1" in seq_record.id:
       referencenames.append("Retro")
    elif "DNA" in seq_record.id:
       referencenames.append("DNA")
    elif ("Helitron" in seq_record.id) or ("helitron" in seq_record.id):
       referencenames.append("DNA")
    elif ("Harbinger" in seq_record.id) or ("harbinger" in seq_record.id):
       referencenames.append("DNA")
    elif ("EnSpm" in seq_record.id):
       referencenames.append("DNA")
    elif ("MuDR" in seq_record.id):
       referencenames.append("DNA")
    elif ("Sola" in seq_record.id) or ("sola" in seq_record.id):
       referencenames.append("DNA")   
    elif ("STOWAWAY" in seq_record.id):
       referencenames.append("DNA")
    elif ("Satellite" in seq_record.id):
       continue
    elif ("Simple_repeat" in seq_record.id):
       continue
    elif ("rRNA" in seq_record.id):
       continue
    else:
       referencenames.append("DNA")
       print "Exceptions:",seq_record.id    
    j=-1
    for nuc in tetranucleotides:
        j+=1
        referencedata[j,i]=dnaseq.count(nuc) # Calculating tetranucleotide frequency
referencedata=np.transpose(referencedata)
print "Reference data is:",referencedata
# Machine learning starts

# Training
clf = svm.SVC(gamma=0.0003, C=1.2)
#clf = svm.SVC(gamma=0.0010, C=2)
clf.fit(querydata, querynames)
referencepredict2=clf.predict(referencedata)
print len(referencenames),len(referencepredict2)
print "Confusion matrix is:",confusion_matrix(referencenames, referencepredict2)
print "Accuracy:", accuracy_score(referencenames, referencepredict2)
#print referencenames,referencepredict
#print "Check\n"
#for i in range(len((seq_records4))):
  #print referencepredict[i],referencenames[i] 

clf = svm.SVC(gamma=0.0010, C=2)
clf.fit(querydata, querynames)
referencepredict=clf.predict(referencedata)
print "Confusion matrix is:",confusion_matrix(referencenames, referencepredict)
print "Accuracy:", accuracy_score(referencenames, referencepredict)

#for i in range(len((seq_records4))):
  #print referencepredict[i],referencenames[i]
