'''
Created on 12 mar. 2021

@author: alba
'''

from Bio import SeqIO
import BCBio
from BCBio import GFF
import sys

def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        if any(fasta) == False:
            print("This file hasn't a fasta format")
            exit()
        #else -> go to the parser function
        else:
            print("***Fasta file***")

def is_gbk(filename):
    with open(filename, "r") as handle:
        genbank = SeqIO.parse(handle, "genbank")
        if any(genbank) == False:
            print("This file hasn't a GenBank format")
        #else -> go to the parser function
        else:
            print("***GenBank file***")

# def is_gff(filename):
#     with open(filename, "r") as handle:

    




if __name__ == "__main__":
    if len(sys.argv) != 2:
        print (__doc__)
        
       
        print ("## Usage format_checker")
        print ("python %s file\n" %sys.argv[0])
       

        sys.exit()

    is_fasta(*sys.argv[1:])
    is_gbk(*sys.argv[1:])
    #is_gff(*sys.argv[1:])
    
    
    