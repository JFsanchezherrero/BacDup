'''
Created on 18 mar. 2021

@author: alba

Modified by JFSanchezHerrero in April 2021
'''

import os
import sys
import re
from Bio import SeqIO
from BCBio import GFF
from builtins import str
from HCGB.functions.aesthetics_functions import debug_message

################################################################################
def is_fasta(filename, debug):
    '''Fasta checker'''
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        
        ## debug messages
        if debug:
            debug_message("Fasta file: ")
            print (filename)
            print (any(fasta))
        
        return any(fasta)

################################################################################
def fasta_seq(filename, debug):
    '''Fasta type checker'''
    seqs = {"dna": re.compile("^[acgt]*$", re.I),
                    "protein": re.compile("^[acdefghiklmnpqrstvwy]*$", re.I)}
    if is_fasta(filename, debug):
        fasta_sequences = SeqIO.parse(open(filename), "fasta")
        for fasta in fasta_sequences:
            sequences = str(fasta.seq)
            dna = seqs["dna"].search(sequences)
            protein = seqs["protein"].search(sequences)
            
            ## debug messages
            if debug:
                debug_message("Type of fasta file: ")
                print ('dna: ' + dna)
                print ('protein: ' + protein)

            return (dna, protein)

################################################################################
def is_gbk(filename, debug):
    '''Genbank format file checker'''
    
    ## debug messages
    if debug:
        debug_message("genbank file?")
                
    for rec in SeqIO.parse(filename, "genbank"):
        if rec.id:
            ## debug messages
            if debug:
                print("True")
            return (True)
        else:
            ## debug messages
            if debug:
                print("False")
            return (False)
                
################################################################################
def is_gff(filename, debug):
    ''' GFF file checker'''
    ## debug messages
    if debug:
        debug_message("GFF file?")

    try:
        with open(filename,"r") as handle:
            limit_info = dict(gff_type=["CDS"])    
            for rec in GFF.parse(handle, limit_info = limit_info):
                if rec.id:
                    ## debug messages
                    if debug:
                        print("True")
                    return (True)
    except:
        ## debug messages
        if debug:
            print("False")
        return (False)
                    
 
################################################################################    
def is_format(filename, debug=False):
    '''Determines input format'''
    
    if debug:
        debug_message("Format for file: " + filename)
    
    ## fasta
    if is_fasta(filename, debug):
        dna, protein = fasta_seq(filename, debug)
        ## DNA
        if dna is not None:
            if debug:
                debug_message("fasta DNA")
            return ("fasta DNA")
        
        ## Protein
        elif protein is not None:
            if debug:
                debug_message("fasta protein" )
            return ("fasta protein")
    
    ## GFF
    if is_gff(filename, debug):
        if debug:
            debug_message("GFF")
        return ("gff")
    
    ## genbank
    if is_gbk(filename, debug):
        if debug:
            debug_message("GenBank")
        return ("gbk")
    
            
    else:
        if debug:
            debug_message("Sorry this file has not any recognizable format:")
            print(os.path.splitext(filename))
        
        return (False)
            

################################################################################
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print (__doc__) 
        print ("## Usage format_checker")
        print ("python %s file\n" %sys.argv[0])
        sys.exit()

    is_format(*sys.argv[1:], debug=True)

    
    