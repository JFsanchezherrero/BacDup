#!/usr/bin/env python3
##############################################################
## Jose F. Sanchez & Alba Moya                              ##
## Copyright (C) 2020-2021                                  ##
##############################################################
'''
Created in March 2021
@author: Jose F. Sanchez-Herrero


This module contains help information for multiple options

'''
from termcolor import colored

import HCGB.sampleParser
import BacDup

##########################
def run_info(options):
    
    ## project help
    if (options.help_project):
        project_help()
        #sampleParser.help_project()
        exit()

    ## help_format option
    if (options.help_format):
        sampleParser.help_format()
        exit()

    ## help_input option
    if (options.input_help):
        BacDup.modules.input_parser.input_help()
        exit()

    
    print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))

##########################
def project_help():
    print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))

##########################
def blast_help():
    print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))
    
    #################################################################
    # qseqid -> query (e.g., unknown gene) sequence id
    # sseqid -> subject (e.g., reference genome) sequence id
    # pident -> percentage of identical matches
    # length -> alignment length (sequence overlap)
    # mismatch -> number of mismatches
    # gapopen -> number of gap openings
    # qstart -> start of alignment in query
    # qend -> end of alignment in query
    # sstart -> start of alignment in subject
    # send -> end of alignment in subject
    # evalue -> expect value
    #            number of expected hits of similar quality (score) that could be found just by chance.
    #            Blast results are sorted by E-value by default (best hit in first line).
    # bitscore -> bit score
    #            requires size of a sequence database in which the current match could be found just by chance.
    #            The higher the bit-score, the better the sequence similarity
    # qlen -> Query sequence length
    # slen -> Subject sequence length
    # aln_pct_qlen -> alignment percentage in query (length/qlen*100)
    # aln_pct_slen -> alignment percentage in subject (length/slen*100)
    ################################################################
    

##########################
def detached_mode():
    print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))

