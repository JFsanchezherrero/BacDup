#!/usr/bin/env python3
##############################################################
## Jose F. Sanchez & Alba Moya                              ##
## Copyright (C) 2020-2021                                  ##
##############################################################
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
from collections import defaultdict 

import HCGB
from HCGB.functions.aesthetics_functions import debug_message
import HCGB.functions.time_functions as time_functions
from HCGB import sampleParser

import BacDup.modules.info as info
import BacDup.scripts.functions as BacDup_functions
import BacDup.scripts.dup_searcher as dup_searcher
import BacDup.modules.input_parser as input_parser

###############
#### TO DO ####
###############
# 6. get results annotation info from annotation csv file ######
#                                          ######################## dup_annot -> dup_searcher -> input_parser  
# 7. classify results by duplicated number #####################
###################################################################################

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
    HCGB.functions.aesthetics_functions.boxymcboxface("Search module")
    print ("--------- Starting Process ---------")
    time_functions.print_time()
    
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
    HCGB.functions.aesthetics_functions.print_sepLine("-",50, False)
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
    start_time_partial = time_functions.timestamp(start_time_total)
    
    ## for each sample
    dict_search_folders = HCGB.functions.files_functions.outdir_project(outdir, arg_dict.project, pd_samples_retrieved, "search", arg_dict.debug)
    annot_search_folders = HCGB.functions.files_functions.outdir_project(outdir, arg_dict.project, pd_samples_retrieved, "dups", arg_dict.debug)

    ## create results
    data2add = pd.DataFrame(columns=("sample", "total_dupGroups", "total_dups", "total_prots"))
    for sample, folder in dict_search_folders.items():
        
        annot_timestamp = os.path.join(annot_search_folders[sample], '.annot_success')
        dup_annot_file = os.path.join(annot_search_folders[sample], 'dup_annot.csv')
        
        ## annotation
        annot_table_file = pd_samples_retrieved.loc[sample, 'annot_table']
            
        if (not HCGB.functions.files_functions.is_non_zero_file(annot_timestamp)):

            ## get results
            file_data = pd_samples_retrieved.loc[sample, 'file_data']
            format = pd_samples_retrieved.loc[sample, 'format']
            filtered_data = dup_searcher.filter_data(sample, file_data, format, arg_dict.pident, arg_dict.evalue, arg_dict.percentage, arg_dict.bitscore, folder, arg_dict.debug)
            
            ## timestamps 
            filter_timestamp = os.path.join(annot_search_folders[sample], '.filter_success')
            if (not HCGB.functions.files_functions.is_non_zero_file(filter_timestamp)):
                #save results as a .csv file
                sort_csv = os.path.abspath(os.path.join(annot_search_folders[sample], 'filtered_results.csv'))
                filtered_data.to_csv(sort_csv, header=True, index=False)
                
                ## print time stamp
                time_functions.print_time_stamp(filter_timestamp)
            else:
                read_time = time_functions.read_time_stamp(filter_timestamp)
                print (colored("\t+ Filter results already available for sample %s [%s]" %(sample, read_time), 'green'))
            
            ## get annotation
            (dup_annot_df, total_dupGroups, total_dups, total_prots) = get_dupannot(filtered_data, annot_table_file, arg_dict.pseudo, arg_dict.debug)
            data2add.loc[sample] = [sample, total_dupGroups, total_dups, total_prots]
            
            ##
            info_dup_file = os.path.join(annot_search_folders[sample], 'info_dup.csv')
            data2add.loc[sample].to_csv(info_dup_file, header=True, index=False)
            
            ## save into file
            dup_annot_df.to_csv(dup_annot_file, header=True)
            
            ## print time stamp
            time_functions.print_time_stamp(annot_timestamp)
            
        else:
            read_time = time_functions.read_time_stamp(annot_timestamp)
            print (colored("\t+ Duplicate annotation already available for sample %s [%s]" %(sample, read_time), 'green'))
    
            ## add info for each
            dup_annot_df = HCGB.functions.main_functions.get_data(dup_annot_file, ',', "index_col=0")
            total_dupGroups = len(dup_annot_df["dup_id"].unique())
            total_dups = dup_annot_df.shape[0]
            annot_table = HCGB.functions.main_functions.get_data(annot_table_file, ',', "index_col=0")
            total_prots = annot_table.shape[0]
            
            ## get data
            data2add.loc[sample] = [sample, total_dupGroups, total_dups, total_prots]
    
    ### report generation
    HCGB.functions.aesthetics_functions.boxymcboxface("Summarizing duplicated search")
    outdir_report = HCGB.functions.files_functions.create_subfolder("report", outdir)

    search_report = HCGB.functions.files_functions.create_subfolder("search", outdir_report)
    
    ## add data2add 
    data2add.to_csv(os.path.join(search_report, 'info_annot.csv'), index=True, header=True)
    
    ## maybe add a summary of the files?
    
    print ("\n*************** Finish *******************")
    start_time_partial = time_functions.timestamp(start_time_total)

    print ("+ Exiting search module.")
    return()

    
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
            HCGB.functions.main_functions.print_all_pandaDF(pd_proteins)
        
            debug_message('pd_annot:', 'yellow')
            HCGB.functions.main_functions.print_all_pandaDF(pd_annot)
        
            debug_message('pd_samples_retrieved:', 'yellow')
            HCGB.functions.main_functions.print_all_pandaDF(pd_samples_retrieved)
        
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
                dict_entries = HCGB.functions.main_functions.file2dictionary(arg_dict.text_file, ',')

                ## check file is readable
                BacDup_functions.file_readable_check(arg_dict.annot_table)
                dict_entries_annot = HCGB.functions.main_functions.file2dictionary(arg_dict.annot_table, ',')

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
                dict_entries = HCGB.functions.main_functions.file2dictionary(arg_dict.fasta_prot, ',')

                ## check file is readable
                BacDup_functions.file_readable_check(arg_dict.annot_table)
                print (colored('\t* Multiple annotation tables provided .......[OK]', 'green'))
                dict_entries_annot = HCGB.functions.main_functions.file2dictionary(arg_dict.annot_table, ',')

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
def get_dupannot(blast_results_df, annot_table_file, pseudo, debug):
    '''Get an annotation information for duplicated proteins'''
    
    ## debug messages
    if debug:
        debug_message('get_dupannot function', 'yellow')
        debug_message('blast_results_df: ', 'yellow')
        print (blast_results_df)
        debug_message('annot_table_file: ' + annot_table_file, 'yellow')
        debug_message('pseudo: ' + str(pseudo), 'yellow')
        
    #get duplicated protein list
    qseqid = list(blast_results_df["qseqid"])
    sseqid =list(blast_results_df["sseqid"])
    prot_id = list(set(qseqid))
    
    ## gets annotation
    annot_table = HCGB.functions.main_functions.get_data(annot_table_file, ',', "index_col=0")
    
    #get filtered_annot table
    filtered_annot = annot_table.loc[prot_id]

    # keep or maintain pseudogenes
    if (pseudo):
        # use them
        print ("+ Pseudogenes would be used in the analysis")
    else:
        filtered_annot = filtered_annot.loc[filtered_annot['pseudo'] != True] 

    ###################################
    ## Get duplicate groups    
    ###################################
    # 1st round
    relations_dict = defaultdict(list) 
    for index, row in blast_results_df.iterrows():
        relations_dict[row['qseqid']].append(row['sseqid'])
    
    ## debug messages
    if debug:
        debug_message('relations_dict: ', 'yellow')
        print (relations_dict)
        
    ## 2nd round
    new_relations_dict = defaultdict(list)
    dups=0
    for key, value in relations_dict.items():
        stop=False
        for dup_id, new_value in new_relations_dict.items():
            if key in new_value:
                stop=True
        if not stop:
            for key2, value2 in relations_dict.items():
                if (key == key2):
                    continue
                else:
                    if (key2 in value): 
                        for i in value2: 
                            if i not in value: 
                                value.extend(i)
            dups += 1
            value.append(key)
            new_relations_dict[str(dups)] = value
            
    ## debug messages
    if debug:
        debug_message('new_relations_dict: ', 'yellow')
        print (new_relations_dict)

    ## Create data
    df_data = pd.DataFrame(columns=('index', 'dup_id'))
    for dup_id, new_value in new_relations_dict.items():
        for i in new_value:
            df_data.loc[i] = (i, dup_id)

    ## merge information    
    dup_annot_df = filtered_annot.join(df_data)
    ## debug messages
    if debug:
        debug_message('dup_annot_df: ', 'yellow')
        HCGB.functions.main_functions.print_all_pandaDF(dup_annot_df)
        
    dup_annot_df = dup_annot_df.drop(columns='index')
    
    ## debug messages
    if debug:
        debug_message('dup_annot_df: ', 'yellow')
        print (dup_annot_df)
        
    ## add info for each
    total_dupGroups = len(dup_annot_df["dup_id"].unique())
    total_dups = dup_annot_df.shape[0]
    total_prots = annot_table.shape[0]
     
    print("#####")
    print("Found %s groups of duplicates with a total of %s proteins duplicated from %s proteins on the original file" % (
        total_dupGroups, total_dups, total_prots))
    print("#####")
    
    return(dup_annot_df, total_dupGroups, total_dups, total_prots)