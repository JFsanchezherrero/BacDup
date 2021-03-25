'''
Created on 18 mar. 2021

@author: alba
'''

import os
from Bio import SeqIO
import BCBio
from BCBio import GFF
import sys


def is_fasta(filename, debug):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        print(any(fasta))
        if any(fasta) == False:
            if(debug):
                print("###")
                print("## DEBUG: This file hasn't a FASTA format ##")
                print("###")
            
        #else -> go to the parser function
        else:
            if(debug):
                print("###")
                print("## DEBUG: FASTA file format ##")
                print("###")
            
        return(any(fasta))
    

def is_gbk(filename, debug):
    with open(filename, "r") as handle:
        genbank = SeqIO.parse(handle, "genbank")
        if any(genbank) == False:
            if(debug):
                print("###")
                print("## DEBUG: This file hasn't a GenBank format ##")
                print("###")
        
        else:
            if(debug):
                print("###")
                print("## DEBUG: GenBank file format ##")
                print("###")
        return (any(genbank))
            
import traceback

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

def is_format (filename, debug=False):
    if is_fasta(filename, debug) == True:
        print("***FASTA file***")
        exit()
    else:             
        if is_gbk(filename, debug) == False:
            print("Sorry this file has not any recognizable format:")
            print(os.path.splitext(filename))
            #Any way to show a preview?
        else:
            print("***GenBank file***")
            

        




if __name__ == "__main__":
    if len(sys.argv) != 2:
        print (__doc__) 
        print ("## Usage format_checker")
        print ("python %s file\n" %sys.argv[0])
        sys.exit()
        
    #is_format(*sys.argv[1:], debug=True)
    is_format(*sys.argv[1:], debug=True)

    