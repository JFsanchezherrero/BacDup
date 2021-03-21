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
import time
from Bio import SeqIO
import HCGB
from HCGB.functions.aesthetics_functions import debug_message
from termcolor import colored

## my modules
import BacDup.scripts.gbf_parser
import BacDup.scripts.gff_parser
import BacDup.scripts.functions as functions

##########################
def input_help():
    print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))
    
    
    compt = {}
    compt["fasta"] = [".fa", ".faa", ".mpfa", ".fna", ".fsa", ".fas", ".fasta"]
    compt["genbank"] = [".genbank", ".gb", ".gbf", ".gbff", ".gbk"]
    compt["GFF"] = [".gff"]
    
    
#############################
def check_annot_file(name, annot_file, Debug):
    """
    This functions checks for each annotation file provided type of input
    and calls appropriate parser: gbf_parser or gff_parser
    """
    ## debug messages
    if (Debug):
        debug_message('+++++++++++++++++++++++++++++++')
        debug_message('check_annot_file function call:', 'yellow')
        debug_message('name: ' + name, 'yellow')
        debug_message('annot_file: ' + annot_file, 'yellow')
        debug_message('+++++++++++++++++++++++++++++++')
    
    ## check file integrity and parse appropriately
    ## file exists non-zero
    ## check format
    ## call parser
    
    #get  annot file absolute path    
    name_file, extension = os.path.splitext(annot_file)
    
    ## GBF
    #prot_file, csv_file = gbf_parser.gbf_parser_caller(arg_dict["annot_file"], output_path, arg_dict["debug"])
    
    ## GFF
    #prot_file, csv_file = gff_parser.gff_parser_caller(arg_dict["annot_file"], arg_dict["ref_file"], output_path, arg_dict["debug"])


##########################
def run_input(arg_dict):
    
    """
    Main function of the input_parser module in BacDup package.
    
    This module prepares data for later gene duplication analysis. 
    
    It allows the user to provide either a single sample, multiple samples or NCBI 
    GenBank IDs to retrieve and obtain the data.    
    """
    
    ## help message
    if (arg_dict.input_help):
        help_input()
        exit()
    
    functions.pipeline_header('BacDup')
    HCGB.functions.aesthetics_functions.boxymcboxface("Preparing input files")
    print ("--------- Starting Process ---------")
    HCGB.functions.time_functions.print_time()
    
    ## init time
    start_time_total = time.time()
    
    ## absolute path for in & out
    #input_dir = os.path.abspath(options.input)
    outdir = os.path.abspath(arg_dict.output_folder)

    ## output folder    
    print ("\n+ Create output folder(s):")
    HCGB.functions.files_functions.create_folder(outdir)

    ## project or detached?
    if arg_dict.detached:
        arg_dict.project = False
        final_dir = outdir
        data_dir = outdir
    else:
        arg_dict.project = True
        print ("+ Generate a directory containing information within the project folder provided")
        final_dir = HCGB.functions.files_functions.create_subfolder("info", outdir)
        data_dir = HCGB.functions.files_functions.create_subfolder("data", outdir)
    
    ## debug messages
    if (arg_dict.debug):
        
        debug_message('+++++++++++++++++++++++++++++++')
        debug_message('Project/Detached option:', 'yellow')
        debug_message('arg_dict.detached: ' + str(arg_dict.detached), 'yellow')
        debug_message('arg_dict.project: ' + str(arg_dict.project), 'yellow')
        debug_message('outdir:' + outdir, 'yellow')
        debug_message('data_dir:' + data_dir, 'yellow')
        debug_message('final_dir:' + final_dir, 'yellow')
        debug_message('+++++++++++++++++++++++++++++++')
        
        
    ## get files
    print ()
    HCGB.functions.aesthetics_functions.print_sepLine("-",50, False)
    print ('+ Getting input information provided... ')
    print ('+ Several options available:')
    print ('\t* Single/Multiple Annotation file:')
    print ('\t  |-- GenBank format files')
    print ('\t  |-- GFF files +  Reference fasta files required')
    print ('\n\t* Single/Multiple NCBI GenBank IDs')
    print ('\n\t* Single/Multiple NCBI taxonomy IDs + Options')
    print ('\n\t* A previous BacDup project folder')
    
    print ('\n+ Check the option provided...')
    time.sleep(1)
    
    #################################################
    ## Parse and obtain the type of input 
    ## information provided
    #################################################
    ## TODO: Now set as mutually_exclusive group. It might be 
    ##       to multiple options
    
    ## --------------------------------------- ##
    ## GFF or GBF file
    ## --------------------------------------- ##
    if (arg_dict.annot_file):
        arg_dict.annot_file = os.path.abspath(arg_dict.annot_file) 
        
        ## multiple files provided
        if (arg_dict.batch):
            if (HCGB.functions.files_functions.is_non_zero_file(arg_dict.annot_file)):
                
                print (colored('\t* Multiple annotation files provided .......[OK]', 'green'))
                dict_entries = HCGB.functions.main_functions.file2dictionary(arg_dict.annot_file, ',')
                
                ## debug messages
                if (arg_dict.debug):
                    debug_message('+++++++++++++++++++++++++++++++')
                    debug_message('Multiple annotation file provided option:', 'yellow')
                    debug_message('dict_entries: ', 'yellow')
                    debug_message(dict_entries, 'yellow')
                    debug_message('+++++++++++++++++++++++++++++++\n\n')
                
                ## for each entry call check_annot_file
                for sample, file_annot in dict_entries.items():
                    check_annot_file(sample, file_annot, arg_dict.debug)
            else:
                print (colored("\n*** ERROR: No readable or accessible file provided. Please check input ***\n", "red"))
                ## debug messages
                if (arg_dict.debug):
                    debug_message('+++++++++++++++++++++++++++++++')
                    debug_message('Multiple annotation file provided option:', 'yellow')
                    debug_message('arg_dict.annot_file: ' + arg_dict.annot_file, 'yellow')
                    debug_message("HCGB.functions.files_functions.is_non_zero_file returned FALSE", 'yellow')
                    debug_message('+++++++++++++++++++++++++++++++\n\n')
                exit()
            
        ### single
        else:
            print (colored('\t* Annotation file:.......[OK]', 'green'))
            check_annot_file(arg_dict.name, arg_dict.annot_file, arg_dict.debug)
            
    ## --------------------------------------- ##
    ## NCBI RefSeq/Genbank IDs: GCA_XXXXXXXX.1; GCF_XXXXXXXXX.1
    ## --------------------------------------- ##
    elif (arg_dict.GenBank_id):
        ## file with multiple
        if (arg_dict.batch):
            arg_dict.GenBank_id = os.path.abspath(arg_dict.GenBank_id)
            
            ## check
            if (HCGB.functions.files_functions.is_non_zero_file(arg_dict.GenBank_id)):
                print (colored('\t* Multiple NCBI GenBank IDs in a file .......[OK]', 'green'))
                
                ## call IDs into a list and create tmp folder
                strains2get = HCGB.functions.main_functions.readList_fromFile(arg_dict.GenBank_id)
                tmp_folder = HCGB.functions.files_functions.create_subfolder("tmp", outdir)
                
                ## debug messages
                if (arg_dict.debug):
                    debug_message('+++++++++++++++++++++++++++++++')
                    debug_message('Multiple NCBI GenBank IDs provided option:', 'yellow')
                    debug_message('arg_dict.GenBank_id: ' + arg_dict.GenBank_id, 'yellow')
                    debug_message('strains2get: ' + strains2get, 'yellow')
                    debug_message('tmp_folder: ' + tmp_folder, 'yellow')
                    debug_message('+++++++++++++++++++++++++++++++')
                
                ## call NCBI_downloader
                data_accID = BacDup.scripts.NCBI_downloader.NCBI_download_list(strains2get, tmp_folder, arg_dict.debug)
                
            else:
                print (colored("\n*** ERROR: No readable or accessible file provided. Please check input ***\n", "red"))
                ## debug messages
                if (arg_dict.debug):
                    debug_message('+++++++++++++++++++++++++++++++')
                    debug_message('GenBank ID file provided option:', 'yellow')
                    debug_message('arg_dict.GenBank_id: ' + arg_dict.GenBank_id, 'yellow')
                    debug_message("HCGB.functions.files_functions.is_non_zero_file returned FALSE", 'yellow')
                    debug_message('+++++++++++++++++++++++++++++++\n\n')
                exit()
        
        ## single GenBank ID
        else:
            print (colored('\t* A NCBI GenBank ID:.......[OK]', 'green'))
            HCGB.functions.aesthetics_functions.print_sepLine("+", 75, False)
            data_accID = NCBIdownload(acc_ID, data_folder, Debug)
    
    ## --------------------------------------- ##
    ## NCBI Taxonomy ID: 
    ## --------------------------------------- ##
    elif (arg_dict.tax_id):
        if (arg_dict.batch):
            print (colored('\t* Multiple NCBI Taxonomy IDs in a file .......[OK]', 'green')) 
        else:
            print (colored('\t* A NCBI Taxonomy ID:.......[OK]', 'green'))

    ## --------------------------------------- ##
    ## Previous BacDup analysis folder
    ## --------------------------------------- ##
    elif (arg_dict.project):
        print (colored('\t* A previous BacDup analysis project folder:.......[OK]', 'green'))
    
    exit()