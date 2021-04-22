#!/usr/bin/env python3

import HCGB
import os
import re
from HCGB.functions.aesthetics_functions import debug_message
from Bio import SeqIO

################################################################################
def columns_rawBLAST_table():
    '''Set the columns to include in a BLAST tabular file'''
    columns = ["qseqid", "sseqid", "pident", "length",
           "mismatch", "gapopen", "qstart", "qend",
           "sstart", "send", "evalue", "bitscore", "qlen", "slen"]
    return(columns)

################################################################################
def columns_accID_table():
    '''Set the columns to include in the accID dataframe'''
    columns= ['new_name','folder','genus','species','taxonomy','genome',
              'annot_file','format_annot_file', 'proteins',
              'plasmids_number','plasmids_ID']
    return(columns)

################################################################################
def columns_annot_table():
    '''Set the columns to include in the annotation tables'''    
    columns = ['rec_id', 'locus_tag', 'protein_id', 'gene', 'start',
               'end', 'strand', 'pseudo', 'product', 'Dbxref', 'inference',
               'EC_number', 'old_locus_tag']
    return(columns)

################################################################################
def file_readable_check(file_given, non_exit=False):
    '''Check file exists and it is ok.
    
    This function checks if file exists and is readable. If not OK by default exits, but if option
    non_exit=True provided, it returns false.
     
    '''
    if (HCGB.functions.files_functions.is_non_zero_file(file_given)):
        return True
    else:
        print (colored("\n*** ERROR: No readable or accessible file provided. Please check input: ***\n " + file_given, "red"))
        if non_exit:
            return (False)
        else:
            exit()
        
################################################################################
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

    
    if (HCGB.functions.files_functions.is_non_zero_file(genome)):
        plasmid_id = []
        for seq_record in SeqIO.parse(genome, "fasta"):
            plasmid_search = re.search(r".*plasmid.*", seq_record.description)
            if plasmid_search:
                name = str( seq_record.id )
                plasmid_id.extend(name)
    else:
        plasmid_id = []
          
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
