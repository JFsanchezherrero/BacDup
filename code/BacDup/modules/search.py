#!/usr/bin/env python3
##############################################################
## Jose F. Sanchez & Alba Moya                              ##
## Copyright (C) 2020-2021                                  ##
##############################################################
from pandas._libs.missing import NA
'''
Created on 13 dic. 2020
@author: alba

Modified in April 2021
@author: Jose F. Sanchez-Herrero

'''

import os
import sys
import argparse
import time
from Bio import SeqIO
from termcolor import colored
import pandas as pd

import HCGB
from HCGB.functions.aesthetics_functions import debug_message
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.main_functions as HCGB_main
from HCGB import sampleParser

import BacDup.modules.info as info
import BacDup.scripts.functions as BacDup_functions
import BacDup.scripts.dup_searcher as dup_searcher
import BacDup.modules.input_parser as input_parser

##########################
def run_search(arg_dict):
    
    """Main function of the search module in BacDup package.
    
    This module searches and create gene duplication analysis. 
    
    It allows the user to provide either a previous parsed data project (NCBI Genbank IDs, taxonomy or user
    annotation data) or a single or multiple samples.    
    """
    
    ## help message
    if (arg_dict.input_help):
        help_input()
        exit()
    
    if (arg_dict.blast_help):
        info.blast_help()
        exit()
    
    if (arg_dict.project_help):
        info.project_help()
        exit()
    
    if (arg_dict.detached_mode_help):
        info.detached_mode()
        exit()
    
    ### Start the analysis
    BacDup_functions.pipeline_header('BacDup')
    HCGB_aes.boxymcboxface("Search module")
    print ("--------- Starting Process ---------")
    HCGB_time.print_time()
    
    ## init time
    start_time_total = time.time()
    
    ## absolute path for in & out
    outdir = os.path.abspath(arg_dict.input_folder)

    ## project or detached?
    if arg_dict.detached:
        arg_dict.project = False
        ## output folder    
        print ("\n+ Create output folder(s):")
        HCGB.functions.files_functions.create_folder(outdir)
    else:
        arg_dict.project = True
    
    ## debug messages
    if (arg_dict.debug):
        debug_message('+++++++++++++++++++++++++++++++')
        debug_message('Project/Detached option:', 'yellow')
        debug_message('arg_dict.detached: ' + str(arg_dict.detached), 'yellow')
        debug_message('arg_dict.project: ' + str(arg_dict.project), 'yellow')
        debug_message('outdir:' + outdir, 'yellow')
        debug_message('+++++++++++++++++++++++++++++++')
        
    ## get files
    print ()
    HCGB_aes.print_sepLine("-",50, False)
    print ('+ Getting information provided... ')
    print ('+ Several options available:')
    print ('\t* BacDup project folder with initiated data')
    print ('\t* Single/Multiple Annotation file:')
    print ('\t  |-- GenBank format files')
    print ('\t  |-- GFF files +  Reference fasta files required')
    print ('\t* Single/Multiple raw BLAST results files')
    print ('\t* Single/Multiple fasta proteins + annotation table')
    
    print ("""\n\n**** NOTE: **** 
    For additional options (e.g. Single/Multiple NCBI GenBank or taxonomy IDs)
    use the input module to accommodate accordingly """)
    time.sleep(1)

    print()
    
    ## parse options
    pd_samples_retrieved = parse_search_options(arg_dict)
   
    ## time stamp
    start_time_partial = HCGB_time.timestamp(start_time_total)
    
    ## for each sample
    dict_search_folders = HCGB.functions.files_functions.outdir_project(outdir, arg_dict.project, 
                                                                        pd_samples_retrieved, "search", arg_dict.debug)
    
    dict_dup_folders = HCGB.functions.files_functions.outdir_project(outdir, arg_dict.project, 
                                                                     pd_samples_retrieved, "dups", arg_dict.debug)

    dict_parse_folders = HCGB.functions.files_functions.outdir_project(outdir, arg_dict.project, 
                                                                     pd_samples_retrieved, "parse", arg_dict.debug)

    ## create results
    data2add = pd.DataFrame(columns=BacDup_functions.columns_dup_table())
    for sample, folder in dict_search_folders.items():
        
        annot_timestamp = os.path.join(dict_dup_folders[sample], '.annot_success')
        dup_annot_file = os.path.join(dict_dup_folders[sample], 'dup_annot.csv')
        
        ## annotation
        annot_table_file = pd_samples_retrieved.loc[sample, 'annot_table']
            
        if (not HCGB.functions.files_functions.is_non_zero_file(annot_timestamp)):

            ## get results
            file_data = pd_samples_retrieved.loc[sample, 'file_data']
            format = pd_samples_retrieved.loc[sample, 'format']
            filtered_data = dup_searcher.filter_data(sample, file_data, format, 
                                                     arg_dict.pident, arg_dict.evalue, 
                                                     arg_dict.percentage, arg_dict.bitscore, 
                                                     folder, arg_dict.debug)
            
            ## timestamps 
            filter_timestamp = os.path.join(dict_dup_folders[sample], '.filter_success')
            if (not HCGB.functions.files_functions.is_non_zero_file(filter_timestamp)):
                #save results as a .csv file
                sort_csv = os.path.abspath(os.path.join(dict_dup_folders[sample], 'filtered_results.csv'))
                filtered_data.to_csv(sort_csv, header=True, index=False)
                
                ## print time stamp
                HCGB_time.print_time_stamp(filter_timestamp)
            else:
                read_time = HCGB_time.read_time_stamp(filter_timestamp)
                print (colored("\t+ Filter results already available for sample %s [%s]" %(sample, read_time), 'green'))
            
            ## get annotation
            (dup_annot_df, data2add_entry) = dup_searcher.get_dupannot(sample, filtered_data, 
                                                                       annot_table_file, arg_dict.debug)
            
            ##
            info_dup_file = os.path.join(dict_dup_folders[sample], 'info_dup.csv')
            data2add_entry.to_csv(info_dup_file, header=True, index=False)
            
            ## save into file
            dup_annot_df.to_csv(dup_annot_file, header=True)
            
            ## print time stamp
            HCGB_time.print_time_stamp(annot_timestamp)
            
        else:
            read_time = HCGB_time.read_time_stamp(annot_timestamp)
            print (colored("\t+ Duplicate annotation already available for sample %s [%s]" %(sample, read_time), 'green'))
    
            ## add info for each
            dup_annot_df = HCGB_main.get_data(dup_annot_file, ',', "index_col=0")
            annot_table = HCGB_main.get_data(annot_table_file, ',', "index_col=0")
            data2add_entry = dup_searcher.get_dup_stats(sample, dup_annot_df, annot_table, arg_dict.debug)
            

        ## add genome length data
        data2add_entry['genome_len'] = ''
        len_df_file = os.path.join(dict_parse_folders[sample], 'length_df.csv')
        if os.path.isfile(len_df_file):
            len_data = HCGB_main.get_data(len_df_file, ',', "header=None")
            data2add_entry['genome_len'] = len_data[1].sum()
        
        ## merge data
        #data2add_entry = data2add_entry.reset_index()
        data2add = data2add.append(data2add_entry, ignore_index=False)
    
    ### report generation
    HCGB_aes.boxymcboxface("Summarizing duplicated search")
    outdir_report = HCGB.functions.files_functions.create_subfolder("report", outdir)
    dups_report = HCGB.functions.files_functions.create_subfolder("dups", outdir_report)
    
    ## add data2add 
    data2add.to_csv(os.path.join(dups_report, 'info_annot.csv'), index=True, header=True)
    
    ## maybe add a summary of the files?
    
    print ("\n*************** Finish *******************")
    start_time_partial = HCGB_time.timestamp(start_time_total)

    print ("+ Exiting search module.")
    return()

##########################
def help_input():
    # TODO: Update this information
    print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))

##########################
def parse_search_options(arg_dict):
    
    ##
    outdir = os.path.abspath(arg_dict.input_folder)

    ## --------------------------------------- ##
    ## Project containing data
    ## --------------------------------------- ##
    if (arg_dict.project):
        print (colored('\t* BacDup project folder:.......[OK]', 'green'))
        
        ## set missing options
        arg_dict.pair=False
        arg_dict.include_all=True
        arg_dict.include_lane=True
        
        ## find samples previously parsed and prepared within a BacDup project structure
        pd_proteins = sampleParser.files.get_files(arg_dict, outdir, "parse", ["fa"], arg_dict.debug)
        pd_proteins = pd_proteins.drop(["dirname", "name", "ext", "tag"], axis=1)
        pd_proteins = pd_proteins.rename(index=str, columns={'sample':'file_data'}) 
        pd_proteins['format'] = 'fasta'
        
        pd_annot = sampleParser.files.get_files(arg_dict, outdir, "parse", ["annot_df.csv"], arg_dict.debug)
        pd_annot = pd_annot.drop(["dirname", "name", "ext", "tag"], axis=1)
        pd_annot = pd_annot.rename(index=str, columns={'sample':'annot_table'}) 
        
        ## merge into pd_samples_retrieved
        pd_samples_retrieved = pd.merge(pd_proteins, pd_annot)
    
        ## debug messages
        if (arg_dict.debug):
            debug_message('pd_proteins:', 'yellow')
            HCGB_main.print_all_pandaDF(pd_proteins)
        
            debug_message('pd_annot:', 'yellow')
            HCGB_main.print_all_pandaDF(pd_annot)
        
            debug_message('pd_samples_retrieved:', 'yellow')
            HCGB_main.print_all_pandaDF(pd_samples_retrieved)
        
    ## --------------------------------------- ##
    ## data on multiple sources
    ## --------------------------------------- ##
    elif (arg_dict.detached):
        print (colored('\t* Detached mode:.......[OK]', 'green'))

        ## parse samples provided
        print()
        
        #########################################################
        ## BLAST raw results provided: either batch or single
        #########################################################
        if (arg_dict.text_file):
            print (colored('\t* BLAST raw results provided:.......[OK]', 'green'))
            print()

            # *************************** ##
            ## Batch file provided
            # *************************** ##
            if (arg_dict.batch):
                ## debug messages
                if (arg_dict.debug):
                    debug_message('+++++++++++++++++++++++++++++++')
                    debug_message('Multiple BLAST results file provided option:', 'yellow')
                    debug_message('arg_dict.text_file: ' + arg_dict.text_file, 'yellow')
    
                ## check if ok
                BacDup_functions.file_readable_check(arg_dict.text_file)
                    
                print (colored('\t* Multiple BLAST results files provided .......[OK]', 'green'))
                dict_entries = HCGB_main.file2dictionary(arg_dict.text_file, ',')

                ## check file is readable
                BacDup_functions.file_readable_check(arg_dict.annot_table)
                dict_entries_annot = HCGB_main.file2dictionary(arg_dict.annot_table, ',')

                ## Check dictionaries contain same information
                if (dict_entries.keys() == dict_entries_annot.keys()):
                    for sample, files in dict_entries.items():
                        ## check annot_table and fasta_file headers are the same ##
                        return_code = dup_searcher.check_annot_table(dict_entries_annot[sample], files, 'BLAST', arg_dict.debug)
                        if not(return_code):
                            print ('Process will continue but sample %s would be discarded' %sample)
                        else:
                            print()
                            ## fill dataframe pd_samples_retrieved
                        
            # *************************** ##
            ## single file provided
            # *************************** ##
            else:
                ## check annot_table and fasta_file headers are the same ##
                return_code = dup_searcher.check_annot_table(arg_dict.annot_table, arg_dict.text_file, 'BLAST', arg_dict.debug)
                if not(return_code):
                    print ('Process will stop here. Please check input files')
                    exit()
                else:
                    print()
                    ## fill dataframe pd_samples_retrieved
                
        #########################################################
        ## annotations file provided: either batch or single
        #########################################################
        elif (arg_dict.annot_file):
            ## debug messages
            if (arg_dict.debug):
                debug_message('Multiple BLAST results file provided option:', 'yellow')
                debug_message('arg_dict.annot_file: ' + arg_dict.annot_file, 'yellow')
            
            ## get input info
            df_accID = input_parser.parse_options(arg_dict)
            if (arg_dict.debug):
                debug_message('df_accID', 'yellow')
                print(df_accID)
            
            ## parse info
            input_parser.parse_information(arg_dict, df_accID, outdir)
            
            ## set missing options
            arg_dict.pair=False
            arg_dict.include_all=True
            arg_dict.include_lane=True
            
            ## find samples previously parsed and prepared within a BacDup project structure
            pd_proteins = sampleParser.files.get_files(arg_dict, outdir, "parse", ["fa"], arg_dict.debug)
            pd_annot = sampleParser.files.get_files(arg_dict, outdir, "parse", ["annot_df.csv"], arg_dict.debug)
        
            ## merge into pd_samples_retrieved
            frames = [pd_proteins, pd_annot]
            pd_samples_retrieved = pd.concat(frames, sort=True, join='outer')
            
            if (arg_dict.debug):
                debug_message('pd_samples_retrieved', 'yellow')
                print(pd_samples_retrieved)
            
        #########################################################
        ## CDS fasta and annotations provided: either batch or single
        #########################################################
        elif arg_dict.fasta_prot:
            
            # *************************** ##
            ## Batch file provided
            # *************************** ##
            if (arg_dict.batch):
                print (colored('\t* Multiple FASTA files provided .......[OK]', 'green'))

                ## debug messages
                if (arg_dict.debug):
                    debug_message('+++++++++++++++++++++++++++++++')
                    debug_message('Multiple Protein FASTA files provided option:', 'yellow')
                    debug_message('arg_dict.fasta_prot: ' + arg_dict.fasta_prot, 'yellow')
    
                ## check if ok
                BacDup_functions.file_readable_check(arg_dict.fasta_prot)
                dict_entries = HCGB_main.file2dictionary(arg_dict.fasta_prot, ',')

                ## check file is readable
                BacDup_functions.file_readable_check(arg_dict.annot_table)
                print (colored('\t* Multiple annotation tables provided .......[OK]', 'green'))
                dict_entries_annot = HCGB_main.file2dictionary(arg_dict.annot_table, ',')

                ## Check dictionaries contain right information
                if (dict_entries.keys() == dict_entries_annot.keys()):
                    for sample, files in dict_entries.items():
                        ## check annot_table and fasta_file headers are the same ##
                        return_code = dup_searcher.check_annot_table(dict_entries_annot[sample], files, 'fasta', arg_dict.debug)
                        if not(return_code):
                            print ('Process will continue but sample %s would be discarded' %sample)
                        else:
                            print()
                            ## fill dataframe pd_samples_retrieved
  
            # *************************** ##
            ## single file provided
            # *************************** ##
            else:
                print (colored('\t* Protein FASTA file provided .......[OK]', 'green'))
                BacDup_functions.file_readable_check(arg_dict.fasta_prot)
                
                ## check file is readable
                print (colored('\t* An annotation table provided .......[OK]', 'green'))
                BacDup_functions.file_readable_check(arg_dict.annot_table)
                
                ## check annot_table and fasta_file headers are the same ##
                return_code = dup_searcher.check_annot_table(arg_dict.annot_table, arg_dict.fasta_prot, 'fasta', arg_dict.debug)
                if not(return_code):
                    print ('Process will stop here. Please check input files')
                    exit()
                else:
                    print()
                    ## fill dataframe pd_samples_retrieved
                    exit()       
   
        ### What??
        else:
            ## Nespresso
            print()
        
    ## return information
    pd_samples_retrieved = pd_samples_retrieved.set_index('new_name')
    return(pd_samples_retrieved)
        
################################################################################

