#!/usr/bin/env python3
##############################################################
## Jose F. Sanchez & Alba Moya                              ##
## Copyright (C) 2020-2021                                  ##
##############################################################
import BacDup
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
import pandas as pd

import BacDup

## my modules
import BacDup.scripts.gbf_parser as gbf_parser
import BacDup.scripts.gff_parser as gff_parser
import BacDup.scripts.format_checker as format_checker
import BacDup.scripts.functions as BacDup_functions
import BacDup.scripts.taxonomy_retrieval as taxonomy_retrieval

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
    
    ## check file integrity and parse appropriately
    ## file exists non-zero
    ## check format
    ## call parser
    if (HCGB.functions.files_functions.is_non_zero_file(annot_file)):
        format = format_checker.is_format(annot_file, Debug)
        
        ## debug messages
        if (Debug):
            debug_message('\nformat_checker.is_format function call:', 'yellow')
            debug_message('format: ' + format, 'yellow')
            
        ## GenBank format
        if (format =='gbk'):
            print()
            ## GBF
            #prot_file, csv_file = gbf_parser.gbf_parser_caller(arg_dict["annot_file"], output_path, arg_dict["debug"])
    
        ## GFF format
        elif (format =='gff'):
            print()
            ## GFF
            #prot_file, csv_file = gff_parser.gff_parser_caller(arg_dict["annot_file"], arg_dict["ref_file"], output_path, arg_dict["debug"])

        ## not valid via this option
        else:
            print()
        

    ## not accessible
    else:
        ## debug messages
        if (Debug):
            debug_message('+++++++++++++++++++++++++++++++')
            debug_message('annot_file:', 'yellow')
            debug_message("HCGB.functions.files_functions.is_non_zero_file returned FALSE", 'yellow')
            
        print (colored("\n*** ERROR: No readable or accessible file provided. Please check input ***\n", "red"))
        exit()

##########################
def run_input(arg_dict):
    
    """Main function of the input_parser module in BacDup package.
    
    This module prepares data for later gene duplication analysis. 
    
    It allows the user to provide either a single sample, multiple samples, NCBI 
    GenBank IDs or NCBI taxonomy IDs to retrieve and obtain the annotation data.    
    """
    
    ## help message
    if (arg_dict.input_help):
        help_input()
        exit()
    
    BacDup_functions.pipeline_header('BacDup')
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

    ## set defaults
    if not (arg_dict.assembly_level):
        arg_dict.assembly_level = 'complete'
    if not (arg_dict.section):
        arg_dict.section = 'genbank'


    ## project or detached?
    if arg_dict.detached:
        arg_dict.project = False
        final_dir = outdir
        data_dir = outdir
    else:
        arg_dict.project = True
        print ("+ Generate a directory containing information within the project folder provided")
        final_dir = HCGB.functions.files_functions.create_subfolder("info", outdir)
    
    ## debug messages
    if (arg_dict.debug):
        
        debug_message('+++++++++++++++++++++++++++++++')
        debug_message('Project/Detached option:', 'yellow')
        debug_message('arg_dict.detached: ' + str(arg_dict.detached), 'yellow')
        debug_message('arg_dict.project: ' + str(arg_dict.project), 'yellow')
        debug_message('outdir:' + outdir, 'yellow')
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
    ## Parse and obtain the type of input information provided
    #################################################
    ## TODO: Now set as mutually_exclusive group. It might be Set to multiple options
    ## ATTENTION: df_accID merge generated dataframe
    
    ## --------------------------------------- ##
    ## GFF or GBF file
    ## --------------------------------------- ##
    if (arg_dict.annot_file):
        arg_dict.annot_file = os.path.abspath(arg_dict.annot_file) 
        
        # *************************** ##
        ## multiple files provided
        # *************************** ##
        if (arg_dict.batch):
            ## debug messages
            if (arg_dict.debug):
                debug_message('+++++++++++++++++++++++++++++++')
                debug_message('Multiple annotation file provided option:', 'yellow')
                debug_message('arg_dict.annot_file: ' + arg_dict.annot_file, 'yellow')

            ## check if ok
            BacDup_functions.file_readable_check(arg_dict.annot_file)
                
            print (colored('\t* Multiple annotation files provided .......[OK]', 'green'))
            dict_entries = HCGB.functions.main_functions.file2dictionary(arg_dict.annot_file, ',')
            
            ## debug messages
            if (arg_dict.debug):
                debug_message('dict_entries: ', 'yellow')
                debug_message(dict_entries, 'yellow')
                debug_message('+++++++++++++++++++++++++++++++\n\n')

        # *************************** ##
        ## single file provided
        # *************************** ##
        else:
            dict_entries = {}
            print (colored('\t* Annotation file:.......[OK]', 'green'))
            if (arg_dict.sample_name):
                sample_name =  arg_dict.sample_name
            else:
                sample_name = "sample"
                
            ##
            dict_entries[sample_name] = arg_dict.annot_file
            
            ## fix dataframe df_accID to match other formats
            ## TODO
        
        ## create dataframe df_accID to match other formats
        df_accID=pd.DataFrame(columns=('new_name','folder','genus','species','taxonomy','genome', 'GFF','GBK', 'proteins','plasmids_number','plasmids_ID'))
        for name, file_annot in dict_entries.items():
            file_annot = os.path.abspath(file_annot)
            
            ## init all
            genome=""
            prot=""
            gff=""
            gbk=""
            plasmid_count = ""
            plasmid_id = "" 
            
            ## debug messages
            if (arg_dict.debug):
                debug_message('+++++++++++++++++++++++++++++++')
                debug_message('dict_entries check annotation files provided option:', 'yellow')
                debug_message('name: ' + name, 'yellow')
                debug_message('file_annot: ' + file_annot, 'yellow')
            
            ## check file is valid
            BacDup_functions.file_readable_check(file_annot)
            
            ## get format
            format = format_checker.is_format(file_annot, arg_dict.debug)
            
            if (arg_dict.debug):
                debug_message('format: ' + format, 'yellow')
            
            ## parse accordingly
            if (format == 'gbk'):
                ## get information from each sample
                (taxonomy, organism) = BacDup.scripts.functions.get_gbk_information(gbk, arg_dict.debug)
                ## plasmid_count, plasmid_id not available
                
            elif (format == 'gff'):
                (taxonomy, organism) = ""
                if (arg_dict.ref_file):
                    arg_dict.ref_file = os.path.abspath(arg_dict.ref_file)
                    BacDup_functions.file_readable_check(arg_dict.ref_file)

                    if (arg_dict.batch):
                        ref_entries = HCGB.functions.main_functions.file2dictionary(arg_dict.ref_file, ',')
                        genome = ref_entries[name]
                    else:
                        genome = arg_dict.ref_file

            ## save into dataframe
            taxonomy_string = ";".join(taxonomy)
            df_accID.loc[len(df_accID)] = (name, dir_path, taxonomy[-1], organism, taxonomy_string, genome, gff, prot, gbk, plasmid_count, ";".join(plasmid_id))

    ## --------------------------------------- ##
    ## NCBI RefSeq/Genbank IDs: GCA_XXXXXXXX.1; GCF_XXXXXXXXX.1
    ## --------------------------------------- ##
    elif (arg_dict.GenBank_id):
        ## get database path
        if (arg_dict.db_folder):
            db_folder = HCGB.functions.files_functions.create_folder(os.path.abspath(arg_dict.db_folder))
        else:
            db_folder = HCGB.functions.files_functions.create_subfolder("db", outdir)

        ## debug messages
        if (arg_dict.debug):
            debug_message('+++++++++++++++++++++++++++++++')
            debug_message('GenBank ID option:', 'yellow')
            debug_message('db_folder: ' + db_folder, 'yellow')
        
        # *************************** ##
        ## batch file
        # *************************** ##
        if (arg_dict.batch):
            arg_dict.GenBank_id = os.path.abspath(arg_dict.GenBank_id)
            
            ## debug messages
            if (arg_dict.debug):
                debug_message('GenBank ID batch file provided:', 'yellow')
                debug_message('arg_dict.GenBank_id: ' + arg_dict.GenBank_id, 'yellow')
            
            ## check is a file and readable
            BacDup_functions.file_readable_check(arg_dict.GenBank_id)
            
            print (colored('\t* Multiple NCBI GenBank IDs in a file .......[OK]', 'green'))
            print()
            
            ## call IDs into a list and create tmp folder
            strains2get = HCGB.functions.main_functions.readList_fromFile(arg_dict.GenBank_id)
                
            ## debug messages
            if (arg_dict.debug):
                debug_message('strains2get: ' + str(strains2get), 'yellow')
                
            # TODO: get info of superkingdom
            group_obtained = 'bacteria'
            
            ## call NCBI_downloader
            df_accID = BacDup.scripts.NCBI_downloader.NCBI_download_list(strains2get, db_folder, arg_dict.debug, 
                                                                         arg_dict.assembly_level, group_given=group_obtained)
            
        # *************************** ##
        ## single GenBank ID
        # *************************** ##
        else:
            ## debug messages
            if (arg_dict.debug):
                debug_message('+++++++++++++++++++++++++++++++')
                debug_message('Single NCBI GenBank IDs provided option:', 'yellow')
                debug_message('arg_dict.GenBank_id: ' + arg_dict.GenBank_id, 'yellow')
                debug_message('db_folder: ' + db_folder, 'yellow')
                debug_message('+++++++++++++++++++++++++++++++')
            
            ## download
            print (colored('\t* A NCBI GenBank ID:.......[OK]', 'green'))
            print()
            HCGB.functions.aesthetics_functions.print_sepLine("+", 75, False)
            df_accID = BacDup.scripts.NCBI_downloader.NCBIdownload(arg_dict.GenBank_id, db_folder, arg_dict.debug)
    
    ## --------------------------------------- ##
    ## NCBI Taxonomy ID: 
    ## --------------------------------------- ##
    elif (arg_dict.tax_id):
        #################
        ## get tax ids
        #################
        if (arg_dict.batch):
            print (colored('\t* Multiple NCBI Taxonomy IDs in a file .......[OK]', 'green'))
            
            ## debug messages
            if (arg_dict.debug):
                debug_message('+++++++++++++++++++++++++++++++')
                debug_message('Multiple NCBI Taxonomy IDs provided option:', 'yellow')
            
            ## check is a file and readable
            BacDup_functions.file_readable_check(arg_dict.tax_id)

            ## get IDs into a list
            taxIDs2get = HCGB.functions.main_functions.readList_fromFile(arg_dict.tax_id)

        else:
            print (colored('\t* A NCBI Taxonomy ID:.......[OK]', 'green'))
            taxIDs2get = [arg_dict.tax_id]
        
        print ()
        
        ##################################
        ## init ete NCBI taxonomy database
        ##################################
        print ('+ Initiate NCBI taxonomy database...')
        ncbi = taxonomy_retrieval.init_db_object(arg_dict.debug)
        
        string_info_total = []
        for taxid in taxIDs2get:
            ## parse
            info = taxonomy_retrieval.parse_taxid(taxid, ncbi, 'unravel', arg_dict.debug)
            print ()
        
            ## debug messages            
            if arg_dict.debug:
                debug_message("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
                debug_message('info\n', "yellow")
                print(info)

            ## append if more        
            string_info_total.extend(info)
            
        ## convert to list of strings
        string_info_total = [str(int) for int in string_info_total]
        
        ## TODO
        ## assuming all belong to same superkingdom
        group_obtained = taxonomy_retrieval.get_superKingdom(string_info_total[0], ncbi, arg_dict.debug)

        #################
        ## get database path
        #################
        if (arg_dict.db_folder):
            db_folder = HCGB.functions.files_functions.create_folder(os.path.abspath(arg_dict.db_folder))
        else:
            db_folder = HCGB.functions.files_functions.create_subfolder("db", outdir)

        ## debug messages            
        if arg_dict.debug:
            debug_message("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
            debug_message('group_obtained: ' + group_obtained, "yellow")
            debug_message('db_folder: ' + db_folder, "yellow")
            debug_message('arg_dict.assembly_level: ' + arg_dict.assembly_level, "yellow")
            debug_message('arg_dict.section: ' + arg_dict.section, "yellow")
            

        ##################################
        ## get GenBank entries selected
        ##################################
        strains2get = taxonomy_retrieval.get_GenBank_ids(db_folder, string_info_total, int(arg_dict.k_random), 
                                                          arg_dict.debug, assembly_level_given=arg_dict.assembly_level,
                                                          group_given=group_obtained, section_given=arg_dict.section)

        ## print list and dictionary of possible and selected taxIDs 
        #HCGB.functions.main_functions.printList2file("./out_file.txt", string_info)
        ## TODO
        #HCGB.functions.main_functions.printDictionary2file("./out_file.txt", string_info)

        #################
        ## call NCBI_downloader
        #################
        df_accID = BacDup.scripts.NCBI_downloader.NCBI_download_list(strains2get, db_folder, arg_dict.debug)

    ## --------------------------------------- ##
    ## Previous BacDup analysis folder
    ## --------------------------------------- ##
    elif (arg_dict.project):
        print (colored('\t* A previous BacDup analysis project folder:.......[OK]', 'green'))
    
        ## create df_accID to store data
    
    
    ### Parse df_accID
    dict_input_folders = HCGB.functions.files_functions.outdir_project(outdir, arg_dict.project, df_accID, "input", arg_dict.debug)

    ## parse each sample retrieved
    for sample, folder_input in dict_input_folders.items():
        ## TODO: Set threads to use in parallel
        print(sample)
        check_annot_file(sample, df_accID, arg_dict.debug)

    
    
    
    
    
    
    exit()