
#Function 1: TRANSLATE Reads in a DNA sequence fasta file,
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
    

#Function 2: Reads in a fasta file of protein sequences,

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

#Function 3: Reads in a fasta file and searches for a set of motifs.

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
