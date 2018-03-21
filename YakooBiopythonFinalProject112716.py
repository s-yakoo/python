
#from Bio import SeqIO

#PART 1 - You will need a Fasta file. One with MULTIPLE SEQUENCE RECORDS.
#   You can use the one in the tutorial (ls_orchid.fasta). 
#   See section 2.4.1  Simple FASTA parsing example

#FUNCTION NAME: ParseFastaFile
#RETURN VALUES: The dictionary with ALL the entries.

#== FUNCTION 1 ==
from Bio import SeqIO
def ParseFastaFile(ls):
    seq={}
    for seq_record in SeqIO.parse(ls, "fasta"):
        list1=seq_record.id.split("|")
        id1=list1[1]
        val=(str(seq_record.seq))
        seq[id1]=val
    return seq

#out=ParseFastaFile("ls_orchid.fasta.txt")
#---------------------------------------------------------------------------------
#PART 2 - You will need a Genbank file. One with MULTIPLE SEQUENCE RECORDS.
#PARAMETERS: 1 (File Name -a string of the name of a GenBank sequence file)
#PURPOSE: The function should:
#           (1) Use Biopython to parse the GenBank file with multiple sequences.
#           (2) Make a dictionary -
#                   (i)  The key will be the record.id
#                   (ii) The value will be a LIST with two items:
#                       1. the description
#                       2. the sequence
#RETURN VALUES: The dictionary with ALL the entries.
#== FUNCTION 2 ==
def ParseGenbankFile(ls2):
    seq={}
    val=[]
    for seq_record in SeqIO.parse(ls2, "genbank"):
        list1=seq_record.id
        val=(seq_record.description,str(seq_record.seq))
        seq[list1]=val
    return seq
    

#---------------------------------------------------------------------------------
#PART 3 -   See sectons 3.8 Transcription AND 3.9  Translation
#PARAMETERS: 1 (A list of sequence nmers)
#PURPOSE: The function should:
#           (1) Iterate through the LIST of DNA sequence nmers. You make them up.
#           (2) Using Biopython, it should transcribe and then translate each sequence.
#           (3) These amino acid sequences should be stored in a list.
#RETURN VALUES: A LIST of all the amino acid sequences (strings)

#== FUNCTION 3 ==
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_rna
def TranscribeTranslate(nmers):
    nmers_list=["AGT","TGG","GTT", "CAT"]
    DNA=""
    protein_list=[]
    for sequence in nmers:
        my_dna=Seq(sequence, generic_dna)
        my_dna=my_dna.transcribe()
        DNA = str(my_dna)
        my_rna=Seq(DNA, generic_rna)
        my_rna=my_rna.translate()
        RNA=str(my_rna)
        protein_list.append(RNA)
    return protein_list
#---------------------------------------------------------------------------------
#PART 4 -  See section 9.6  EFetch: Downloading full records from Entrez
#       This function will use a pubmed ID number to fetch a file from the
#       internet
#       and save it to a file on your computer, under the pubmed id name.

#FUNCTION NAME: FetchGenbankFile
#PARAMETERS: 1 (A string - Pubmed ID number)
#PURPOSE: The function should:
#           (1) Get the genbank file as text from the nucleotide database.
#           (2) Write it to a file. (See below)
#RETURN VALUES: 1 (True)

#== FUNCTION 4 ==
from Bio import Entrez
def FetchGenbankFile(id_Number):
    file_name=""
    Entrez.email = "syakoo@sdsu.edu"
    handle = Entrez.efetch(db="nucleotide", id=id_Number, rettype="gb", retmode="text")
    file_name = str(id_Number) + ".gbk.txt"
    fout=open(file_name,'w')
    fout.write(handle.read())
    fout.close()
    return 1
