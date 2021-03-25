#!/usr/bin/env python3

import HCGB
import os
import re
from HCGB.functions.aesthetics_functions import debug_message
from Bio import SeqIO


##########################################################################################
def get_gbk_information(gbk, debug):
    ## read Genbank file to retrieve information for each samle
    ## https://biopython.org/wiki/SeqRecord
    # get
    for index, record in enumerate(SeqIO.parse(gbk, "genbank")):
        
        if (index == 0): ## only for first entry == Main chromosome
            if debug:
                debug_message("******************************************")
                debug_message("SeqIO.read(gbk, 'genbank') info:", color="yellow")
                debug_message("record", color="yellow")
                debug_message(record, color="yellow")
                
            organism = record.annotations['source']
            taxonomy = record.annotations['taxonomy']
    
    ##
    return(taxonomy, organism)

##########################################################################################
def get_files_annotation(folder, debug):
    '''
    Code retrieve from BacterialTyper database_generator.py script
    '''
    ## check if files are gunzip
    files = os.listdir(folder)
    genome=""
    prot=""
    gff=""
    gbk=""
    for f in files:
        if f.endswith('.fna'):
            genome = os.path.join(folder, f)
        elif f.endswith('.gff'):
            gff = os.path.join(folder, f)
        elif f.endswith('.gbk'):
            gbk = os.path.join(folder, f)
        elif f.endswith('.gbff'):
            gbk = os.path.join(folder, f)
        elif f.endswith('.faa'):
            prot = os.path.join(folder, f)

    ## debug messages
    if debug:
        debug_message("-----------------------------------------")
        debug_message("Return info get_files_download", color="yellow")
        debug_message("genome: " + genome, color="yellow")
        debug_message("prot: " + prot, color="yellow")
        debug_message("gff: " + gff, color="yellow")
        debug_message("gbk: " + gbk, color="yellow")

    return(genome, prot, gff, gbk)

########################################
def get_plasmids(genome, Debug):
    """
    Parses genome fasta file and retrieves plasmid sequences IDs        
    """

    plasmid_id = []
    for seq_record in SeqIO.parse(genome, "fasta"):
        plasmid_search = re.search(r".*plasmid.*", seq_record.description)
        if plasmid_search:
            name = str( seq_record.id )
            plasmid_id.append(name)
    
    ##
    return(len(plasmid_id), plasmid_id)
    
########################################
def pipeline_header(option):
    """
    Prints a common header for the pipeline including name, author, copyright and year.        
    """
    
    print ("\n")
    HCGB.functions.aesthetics_functions.print_sepLine("#", 70, False)
    
    print('#', '{: ^66}'.format("BacDup pipeline"), '#')
    print('#', '{: ^66}'.format("Jose F. Sanchez & Alba Moya Garces"), '#')
    print('#', '{: ^66}'.format("Copyright (C) 2020-2021"), '#')
    HCGB.functions.aesthetics_functions.print_sepLine("#", 70, False)
