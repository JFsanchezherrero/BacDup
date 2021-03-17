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
            
        #else -> go to the parser function
        else:
            print("***Fasta file***")
            
        return (any(fasta))

def is_gbk(filename):
    with open(filename, "r") as handle:
        genbank = SeqIO.parse(handle, "genbank")
        if any(genbank) == False:
            print("This file hasn't a GenBank format")
        
        #else -> go to the parser function
        else:
            print("***GenBank file***")
        return (any(genbank))
            

def is_gff(filename):
    try:
        with open(filename,"r") as handle:
            GFF.parse(handle)
            print("***GFF file***")
            return(True)
    except:
        print("This file hasn't a GFF format")
        return(False)
#     with open(filename, "r") as handle:

    




if __name__ == "__main__":
    if len(sys.argv) != 2:
        print (__doc__)
        
       
        print ("## Usage format_checker")
        print ("python %s file\n" %sys.argv[0])
       

        sys.exit()

    is_fasta(*sys.argv[1:])
    is_gbk(*sys.argv[1:])
    is_gff(*sys.argv[1:])
    #is_gff(*sys.argv[1:])
    
    
    