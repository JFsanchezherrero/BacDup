#!/usr/bin/env python3
##############################################################
## Jose F. Sanchez & Alba Moya                              ##
## Copyright (C) 2020-2021                                  ##
##############################################################
from pickle import FALSE
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

import BacDup.scripts.blast_caller as blast_caller
import BacDup.scripts.functions as BacDup_functions
from HCGB.functions.aesthetics_functions import debug_message

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
def get_data(file_data, format, out_folder, debug):
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
        
        ## FIXME
        raw_blast = blast_caller.create_blast_results(arg_dict)
        raw_blast = pd.read_csv(raw_blast, sep="\t", header = None, names=columns)
        
    return (raw_blast)

################################################################################
def filter_data(file_data, format, pident, evalue, percent, bitscore, out_folder, debug):
    
    '''
    Get duplicated proteins and generate a filtered and sort dataframe depending on arguments values
    if a BLAST results text file is provided, it should have next columns names.    
    '''
    ## get BLAST results and filter
    raw_blast = get_data(file_data, format, out_folder, debug)
    filtered_results = filter_BLAST(raw_blast, pident, evalue, percent, bitscore)

    #sort by aln_pct_qlen (desc), evalue(asc), bitscore(desc)
    by_alnpct = filtered_results.sort_values(["aln_pct_qlen", "evalue", "bitscore"],
                                   ascending=[False, True, False])

    #save results as a .csv file
    sort_csv_file = os.path.abspath(os.path.join(out_folder, 'filtered_results.csv'))
    by_alnpct.to_csv(sort_csv, header=True, index=False)
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
