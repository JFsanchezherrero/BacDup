'''
Created on 24 mar. 2021

@author: alba
'''
import os
from Bio import SeqIO
import BCBio
from BCBio import GFF
import sys

#filename = "/home/alba/git/BacDup/developer/test_data/example_denovo/example_genomic.fna"

def is_fasta(filename, debug):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)

def is_gbk(filename, debug):
    with open(filename, "r") as handle:
        genbank = SeqIO.parse(handle, "genbank")
        return any(genbank)
    
def is_format (filename, debug=False):
    if is_fasta(filename, debug) == True:
        print("fasta")
        exit()
    else:             
        if is_gbk(filename, debug) == False:
            print("Sorry this file has not any recognizable format:")
            print(os.path.splitext(filename))
            #Any way to show a preview?
        else:
            print("gbk")
            
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print (__doc__) 
        print ("## Usage format_checker")
        print ("python %s file\n" %sys.argv[0])
        sys.exit()
        
    #is_format(*sys.argv[1:], debug=True)
    is_format(*sys.argv[1:], debug=True)