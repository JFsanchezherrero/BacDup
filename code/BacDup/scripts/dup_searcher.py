#!/usr/bin/env python3
##############################################################
## Jose F. Sanchez & Alba Moya                              ##
## Copyright (C) 2020-2021                                  ##
##############################################################
'''
Created on 3 dic. 2020
@author: alba

Modified in April 2021
@author: Jose F. Sanchez-Herrero
'''
import os
import sys
import pandas as pd
import numpy as np
from termcolor import colored
from Bio import SeqIO, Seq

import argparse
from argparse import ArgumentParser

import HCGB
from HCGB.functions.aesthetics_functions import debug_message
import HCGB.functions.time_functions as time_functions

import BacDup
import BacDup.scripts.functions as BacDup_functions

################################################################################
def create_blast_results(sample, fasta_file, outdir, debug):
    '''Creates BLAST results for each fasta vs. itself'''
    
    #phr is the header file, pin is the index file, psq is the sequence file
    
    ## debug messages
    if debug:
        debug_message('create_blast_results function call:', 'yellow')
        debug_message('sample: ' + sample, 'yellow')
        debug_message('fasta_file: ' + fasta_file, 'yellow')
        debug_message('outdir: ' + outdir, 'yellow')
    
    ## output file
    raw_blast = os.path.abspath(os.path.join(outdir, "BLAST_raw_results.tsv"))

    ## timestamps 
    db_timestamp = os.path.join(outdir, '.db_success')
    search_timestamp = os.path.join(outdir, '.blast_success')
        
    if (not HCGB.functions.files_functions.is_non_zero_file(db_timestamp) and not HCGB.functions.files_functions.is_non_zero_file(search_timestamp)):

        ## get binaries
        (makeblastdb_exe, blastp_exe) = BacDup.modules.config.get_exe('BLAST', debug)
        makeblastdb_exe = "/usr/bin/makeblastdb" 
        blastp_exe = "/usr/bin/blastp"
        
        ## check if db is indexed already
        db_path_name = os.path.join(os.path.abspath(outdir), sample + '_db')
        if (not HCGB.functions.files_functions.is_non_zero_file(db_timestamp)):
            ## generate blastdb for genome
            HCGB.functions.blast_functions.makeblastdb(db_path_name, fasta_file, makeblastdb_exe, 'prot') # HCGB function    
        
            ## print time stamp
            time_functions.print_time_stamp(db_timestamp)
        
        else:
            print (colored("\t+ BLAST database already available for sample %s [%s]" %(sample, read_time), 'green'))
            
        ## create blastp outfile
        HCGB.functions.blast_functions.blastp(blastp_exe, raw_blast, db_path_name, fasta_file, 1) # HCGB function

        ## print time stamp
        time_functions.print_time_stamp(search_timestamp)
    else:
        read_time = time_functions.read_time_stamp(search_timestamp)
        print (colored("\t+ Duplicate search already available for sample %s [%s]" %(sample, read_time), 'green'))
            
    return (raw_blast)

################################################################################
def check_annot_table(annot_table, file, format, debug):
    '''Check annotation table provided matches the BLAST or protein fasta file provided'''
    
    ## read annotation
    annotation_table = pd.read_csv(annot_table, sep=",", index_col=0 )
    #annotation_table = annotation_table.drop("Unnamed: 0",axis=1)
    
    ## debug messages
    if debug:
        debug_message('check_annot_table function:', 'yellow')
        debug_message("list(annotation_table.columns)", 'yellow')
        print(list(annotation_table.columns))
        debug_message("BacDup_functions.columns_annot_table()", 'yellow')
        print(BacDup_functions.columns_annot_table())
    
    ## skip if not OK
    if not (list(annotation_table.columns) == BacDup_functions.columns_annot_table()):
        print('ERROR: Annotation table does NOT match desired input format')
        return(False)
    
    ## read BLAST sequence results
    if format=="blast":
        ## TODO
        print()
    
    ## protein fasta file
    elif(format=="fasta"):
        with open (file) as in_handle:
            ref_recs = SeqIO.to_dict(SeqIO.parse(in_handle, "fasta"))

        protein_IDs = list(ref_recs.keys()) ## get sequence headers
        index_annot = list(annotation_table.index) ## get annotation labels
        
        ## debug messages
        if debug:
            debug_message("list(annotation_table.index).sort()", 'yellow')
            print(index_annot)
            debug_message("list(ref_recs.keys()).sort()", 'yellow')
            print(protein_IDs)
        
        if (index_annot == protein_IDs):
            return(True)
        else:
            return(False)
        
################################################################################
def filter_BLAST(raw_blast, pident, evalue, bitscore, percent):
    '''Filter BLAST tabular results
    Created on 3 dic. 2020
    @author: alba
    '''    
    
    #fill aln_pct columns
    raw_blast["aln_pct_qlen"] = (raw_blast["length"]/raw_blast["qlen"]*100).round(2)
    raw_blast["aln_pct_slen"] = (raw_blast["length"]/raw_blast["slen"]*100).round(2)

    #filter mirror pairs
    del_mirror =pd.DataFrame(np.sort(raw_blast[["qseqid", "sseqid"]], axis=1))
    raw_blast = raw_blast[~del_mirror.duplicated()]
          
    ## filter using different parameters
    filter_pident= raw_blast["pident"] >= pident
    filter_evalue = raw_blast["evalue"] <= evalue
    filter_bitscore = raw_blast["bitscore"] >= bitscore
    filter_alnpct_slen = raw_blast["aln_pct_slen"] >= percent
    filter_alnpct_qlen = raw_blast["aln_pct_qlen"] >= percent
    filter_id = raw_blast["qseqid"] != raw_blast["sseqid"]

    filtered_results = raw_blast[filter_pident & filter_evalue & filter_bitscore & filter_alnpct_qlen & filter_alnpct_slen & filter_id]
    
    return(filtered_results)

################################################################################
def get_data(sample, file_data, format, out_folder, debug):
    '''Function to get BLAST results'''
    
    file_data = os.path.abspath(file_data)
        
    ## debug messages
    if (debug):
        debug_message('dup_searcher.get_data:', 'yellow')
        debug_message('file_data:' + file_data, 'yellow')
        debug_message('format:' + format, 'yellow')
    
    ## check file is readable
    BacDup_functions.file_readable_check(file_data)
    
    ## parse accordingly
    if (format=='blast_raw'):
        ## FIXME
        raw_blast = pd.read_csv(file_data, sep="\t", header = None, names=BacDup_functions.columns_rawBLAST_table())
                
    elif (format=='fasta'):
        
        raw_blast = create_blast_results(sample, file_data, out_folder, debug)
        raw_blast = pd.read_csv(raw_blast, sep="\t", header = None, names=BacDup_functions.columns_rawBLAST_table())
        
    return (raw_blast)

################################################################################
def filter_data(sample, file_data, format, pident, evalue, percent, bitscore, out_folder, debug):
    
    '''
    Get duplicated proteins and generate a filtered and sort dataframe depending on arguments values
    if a BLAST results text file is provided, it should have next columns names.    
    '''
    ## get BLAST results and filter
    raw_blast = get_data(sample, file_data, format, out_folder, debug)
    filtered_results = filter_BLAST(raw_blast, pident, evalue, percent, bitscore)

    #sort by aln_pct_qlen (desc), evalue(asc), bitscore(desc)
    by_alnpct = filtered_results.sort_values(["aln_pct_qlen", "evalue", "bitscore"],
                                       ascending=[False, True, False])
    return(by_alnpct)

###############################################################
def help_options():
    print ("\nUSAGE: python %s file format pident--evalue--percent--bitscore Debug[True/False]\n"  %os.path.realpath(__file__))

############################################################
def main():
    ## this code runs when call as a single script

    ## control if options provided or help
    if len(sys.argv) > 1:
        print ("")
    else:
        help_options()
        exit()
        
    file_data = os.path.abspath(sys.argv[1])
    format = sys.argv[2]
    filters = sys.argv[3].split('--')
    debug = sys.argv[4]
        
    filter_data(file_data, format, pident, evalue, percent, bitscore, '.', debug)

############################################################
if __name__== "__main__":
    main()
