#!/usr/bin/env python3
##############################################################
## Jose F. Sanchez & Alba Moya                              ##
## Copyright (C) 2020-2021                                  ##
##############################################################
'''
Created on 25 oct. 2020
@author: alba

Modified in March 2021
@author: Jose F. Sanchez-Herrero
'''

############
    # 1: identificar annot_file: GFF / Genbank / None
    
    # 2: Si es genbank: gbf_parser
    
    # 3: Si es GFF:
    # 3.1: Comprobar ref_file: Yes/No
    # 3.2: gff_parser(annot_file, ref_file)
    
    # ToDO: create gtf parser
#############


## useful imports
import os
import sys
import argparse
from Bio import SeqIO
import HCGB
from termcolor import colored

## my modules
import BacDup.scripts.gbf_parser
import BacDup.scripts.gff_parser

##########################
def help_input():
    print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))

##########################
def run_input(arg_dict):
    
    """
    Main function of the input_parser module in BacDup package.
    
    This module prepares data for later gene duplication analysis. 
    
    It allows the user to provide either a single sample, multiple samples or NCBI 
    GenBank IDs to retrieve and obtain the data.    
    """
    
    ## help message
    if arg_dict.help_input():
        help_input()
        exit()
    
    compt = {}
    compt["fasta"] = [".fa", ".faa", ".mpfa", ".fna", ".fsa", ".fas", ".fasta"]
    compt["genbank"] = [".genbank", ".gb", ".gbf", ".gbff", ".gbk"]
    compt["GFF"] = [".gff"]
    
    if arg_dict["debug"]:
        print("## DEBUG: Format input file ")
        print(compt)
        
    ## output
    output_path= os.path.abspath(arg_dict["out_folder"]) 
    HCGB.functions.file_functions.create_folder(output_path)
     
    #get  annot file absolute path    
    file_name_abs_path = os.path.abspath(arg_dict["annot_file"])
    name_file, extension = os.path.splitext(file_name_abs_path)
    
    if arg_dict["debug"]:
        print("## DEBUG: name_file and extension ")
        print(os.path.splitext(file_name_abs_path))
         
    #get protein fasta file and annotation csv file    
    if extension in compt["genbank"]:
        if not arg_dict["ref_file"] is None:
            print("######")
            print("Sorry, we don't need a ref file to parse a genbank file")
            print("######")
        else:
        ## call gbf_parser
            prot_file, csv_file = gbf_parser.gbf_parser_caller(arg_dict["annot_file"], output_path, arg_dict["debug"])
            return (prot_file, csv_file)
        if arg_dict["debug"]:
                print("## DEBUG: prot_file ")
                print(prot_file)
    
    elif extension in compt["GFF"]:
        if arg_dict["ref_file"]==None:
            print("######")
            print("Please provide a ref_file FASTA format")
            print("######")
            print(parser.print_help())
        else:
            ref_name_abs_path = os.path.abspath(arg_dict["ref_file"])
            ref_name, ref_extension = os.path.splitext(ref_name_abs_path)
            
            if arg_dict["debug"]:
                print("## DEBUG: ref_name and extension ")
                print(os.path.splitext(ref_name_abs_path))
            
            if ref_extension in compt["fasta"]:
                #protein_gff.protein_recs(arg_dict["annot_file"], arg_dict["ref_file"])
                prot_file, csv_file = gff_parser.gff_parser_caller(arg_dict["annot_file"], arg_dict["ref_file"], output_path, arg_dict["debug"])
                return (prot_file, csv_file)
    
            else:
                print("######")
                print("Compatible ref_file formats: ")
                print(compt["fasta"])
                print("######")
    else:
        print("######")
        print("Compatible annot_file formats: ")
        print("#GenBank: ")
        print(compt["genbank"])
        print("#GFF3: ")
        print(compt["GFF"])
        print("######")