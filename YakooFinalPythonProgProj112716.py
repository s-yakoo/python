##NAME: Sylvia Yakoo

#FINAL PROJECT

##SPECIAL NOTES: DO NOT INCLUDE ANY FUNCTION CALLS IN YOUR FINAL CODE
## ALSO, USE THE FUNCTION NAMES EXACTLY AS I HAVE WRITTEN THEM
##  AND DO NOT CHANGE THE NUMBER OR TYPES OF PARAMETERS

#Function 1: TRANSLATE Reads in a DNA sequence fasta file,
#and outputs a Fasta file of protein translations and it
#should be able to do standard (eukaryotic) translation only.
# In other words, don't use the mitochondrial code and assume that the
# DNA you have is eukaryotic protein coding DNA.

#(You don't have to change the names of the fasta sequence titles.)
standard_code = {
     "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S",
     "UCC": "S", "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y",
     "UAA": "*", "UAG": "*", "UGA": "*", "UGU": "C", "UGC": "C",
     "UGG": "W", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
     "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAU": "H",
     "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R", "CGC": "R",
     "CGA": "R", "CGG": "R", "AUU": "I", "AUC": "I", "AUA": "I",
     "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
     "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGU": "S",
     "AGC": "S", "AGA": "R", "AGG": "R", "GUU": "V", "GUC": "V",
     "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A",
     "GCG": "A", "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
     "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

def DNA2Prot(f1, f2="translated_fasta.txt"):
    fn=open(f1,'U')
    fout=open(f2,'w')
    for line in fn:
        if ">" in line:
            fout.write(line)
            pass
        else:
            codons=[]
            i=0
            protein=""
            line=line.replace("T","U")
            line=line.strip()
            for i in range(0,len(line),3):
                codon=line[i:i+3]
                if len(codon)== 3:
                    codons.append(codon)
                    amino_acid = standard_code[codon]
                    protein = protein + amino_acid
                else:
                    pass
            fout.write(protein)
            fout.write("\n")

    fn.close()
    fout.close()
    return 1
    
#out=DNA2Prot("testDNAseq2016.txt","translated_fasta.txt")

#Function 2: Reads in a fasta file of protein sequences,
#makes a table of amino acids frequencies for each protein
#sequence and outputs it to an excel readable file (tab delimited).


def AAfreq(f1,f2="aatable.xls"):
    fn = open(f1, 'r')
    fout = open(f2, 'w')
    count = 0
    countA = 0
    fseq = fn.readlines()
    fdict = {}
    for line in fseq:
        if ">" in line:
            fseq.remove(line)
    for line in fseq:
        for line2 in fseq(countA):
            if (line2 != "\n"):
                if line2 not in fdict:
                    fdict[line2] = 1
                else:
                    fdict[line2] += 1
        countA += 1
        if (count == 0):
            fout.write("AminoAcid\t")
            for line in fseq:
                fout.write(str(line) + "\t")
            fout.write("\n")
        fout.write("Seq" +str(count) + "\t")
        fout.write("\n")
        count += 1
    fout.close()

    return f2

#out=AAfreq("testProtein.txt","aatable.xls")



#Function 3: Reads in a fasta file and searches for a set of motifs.
# The function should read in a LIST of motifs and outputs a
#tab-delimited file of the number of times that each motif was found
#in each sequence. The file should open nicely in Excel or a spreadsheet program.



import re

def MotifFinder(f1, motif, f2="motifs.xls"):
    fn = open(f1, 'r')
    fout = open(f2, 'w')
    count = 1
    fseq = fn.readlines()
    for line in fseq:
        if ">" in line:
            fseq.remove(line)
    fout.write("SeqName" + "\t" + "M1" + "\t" + "Hits" + "\t" + "M2" + "\t" + "Hits""\n")
    for line in fseq:
        fout.write("Seq" +str(count))
        count += 1
        for line2 in motif:
            motifs=re.findall(line2,line)
            fout.write("\t")
            fout.write(line2)
            fout.write("\t")
            fout.write(str(len(motifs)))
        fout.write("\n")
    fout.close()
    return f2

        
            
#Example Function Call:
#MotifFinder("testProtein.txt",["MN[A-Z]","V[A-Z]R[ML]"])
