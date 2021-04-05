'''
Created on 18 mar. 2021

@author: alba
'''

import os
from Bio import SeqIO
import BCBio
from BCBio import GFF
import sys
from builtins import str
import re




def is_fasta(filename, debug):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)

def fasta_seq(filename, debug):
    seqs = {"dna": re.compile("^[acgt]*$", re.I),
                    "protein": re.compile("^[acdefghiklmnpqrstvwy]*$", re.I)}
    if is_fasta(filename, debug):
        fasta_sequences = SeqIO.parse(open(filename), "fasta")
        for fasta in fasta_sequences:
            sequences = str(fasta.seq)
            dna = seqs["dna"].search(sequences)
            protein = seqs["protein"].search(sequences)
            return (dna, protein)

def is_gbk(filename, debug):
    with open(filename, "r") as handle:
        genbank = SeqIO.parse(handle, "genbank")
        return any(genbank)
            


## FIXME!!

def is_gff(filename):
    with open(filename,"r") as handle:
        try:
            GFF.parse(handle)
            print("***GFF3 file***")
            return(True)
        except AssertionError:
            print("This file hasn't a GFF format")
            return(False)
#     with open(filename, "r") as handle:

#### TODO import hbcg debug
## TODO comprobar columnas csv en el script "functions"

#filename = "/home/alba/git/BacDup/developer/test_data/example_denovo/example_genomic.fna"


    
def is_format (filename, debug=False):
    if is_fasta(filename, debug):
        dna, protein = fasta_seq(filename, debug)
        if dna is not None:
            if debug:
                print("## DEBUG: fasta DNA ##")
            return str("fasta DNA")
        elif protein is not None:
            if debug:
                print("## DEBUG: fasta protein ##")
            return str("fasta protein")
        exit()
    elif is_gbk(filename, debug):
        if debug:
            print("## DEBUG: gbk ##")
        return str("fasta")
        exit()
        
    else:
        print("## DEBUG: Sorry this file has not any recognizable format:")
        print(os.path.splitext(filename))
        print("###")
        #Any way to show a preview?
            

            
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print (__doc__) 
        print ("## Usage format_checker")
        print ("python %s file\n" %sys.argv[0])
        sys.exit()

    is_format(*sys.argv[1:], debug=True)

    
    