#!/usr/bin/env python3
##############################################################
## Jose F. Sanchez & Alba Moya                              ##
## Copyright (C) 2020-2021                                  ##
##############################################################
'''
Created on 8 june 2021
@author: Jose F. Sanchez-Herrero
'''

import os
import sys
import argparse
import time
#from Bio import SeqIO
from termcolor import colored
import pandas as pd

import HCGB
from HCGB.functions.aesthetics_functions import debug_message
import HCGB.functions.time_functions as time_functions
from HCGB import sampleParser

import BacDup.modules.info as info
import BacDup.scripts.functions as BacDup_functions

##########################
def run_report(arg_dict):
    
    """Main function of the plot and report generation module in BacDup package.
    
    This module searches and create gene duplication analysis. 
    
    It allows the user to provide either a previous parsed data project (NCBI Genbank IDs, taxonomy or user
    annotation data) or a single or multiple samples.    
    """
    
    ## help message
    if (arg_dict.input_help):
        help_input()
        exit()
    
    if (arg_dict.project_help):
        info.project_help()
        exit()
    
    ### Start the analysis
    BacDup_functions.pipeline_header('BacDup')
    HCGB.functions.aesthetics_functions.boxymcboxface("Report generation module")
    print ("--------- Starting Process ---------")
    time_functions.print_time()
    
    ## init time
    start_time_total = time.time()
    
    ## absolute path for in & out
    outdir = os.path.abspath(arg_dict.input_folder)

    ## default
    arg_dict.project = True
    arg_dict.batch = False
    
    ## debug messages
    if (arg_dict.debug):
        debug_message('+++++++++++++++++++++++++++++++')
        debug_message('Project/Detached option:', 'yellow')
        debug_message('arg_dict.detached: Not available', 'yellow')
        debug_message('arg_dict.project: ' + str(arg_dict.project), 'yellow')
        debug_message('outdir:' + outdir, 'yellow')
        debug_message('+++++++++++++++++++++++++++++++')

    ## get files
    print ()
    HCGB.functions.aesthetics_functions.print_sepLine("-",50, False)
    print ('+ Getting information provided... ')
    time.sleep(1)
    print()
    
    ## parse options
    pd_samples_dups = sampleParser.files.get_files(arg_dict, outdir, "dups", ["dup_annot.csv"], arg_dict.debug)
    pd_samples_dups = pd_samples_dups.drop(["dirname", "name", "ext", "tag"], axis=1)
    pd_samples_dups = pd_samples_dups.rename(index=str, columns={'sample':'file_data'}) 
    pd_samples_dups['format'] = 'dup_annot'

    ## debug messages
    if (arg_dict.debug):
        debug_message("pd_samples_dups",'yellow')
        HCGB.functions.main_functions.print_all_pandaDF(pd_samples_dups)
        
    pd_info = sampleParser.files.get_files(arg_dict, outdir, "parse", ["length_df.csv"], arg_dict.debug)
    pd_info = pd_info.drop(["dirname", "name", "ext", "tag"], axis=1)
    pd_info = pd_info.rename(index=str, columns={'sample':'length_table'}) 
    
    ## debug messages
    if (arg_dict.debug):
        debug_message("pd_info",'yellow')
        HCGB.functions.main_functions.print_all_pandaDF(pd_info)
    
    ## merge into pd_samples_retrieved
    pd_samples_retrieved = pd.merge(pd_samples_dups, pd_info)

    ## debugging messages
    if arg_dict.debug:
        debug_message("pd_samples_retrieved", 'yellow')
        HCGB.functions.main_functions.print_all_pandaDF(pd_samples_retrieved)
    
    ## time stamp
    start_time_partial = time_functions.timestamp(start_time_total)
    
    ## for each sample
    dict_plot_folders = HCGB.functions.files_functions.outdir_project(outdir, arg_dict.project, 
                                                                      pd_samples_retrieved, "dup_plot", arg_dict.debug)
    
    ## debugging messages
    if arg_dict.debug:
        debug_message("dict_plot_folders", 'yellow')
        print (dict_plot_folders)
    
    ## create results
    for sample, folder in dict_plot_folders.items():
        plot_timestamp = os.path.join(dict_plot_folders[sample], '.plot_success')
        dup_annot_file = pd_samples_retrieved[pd_samples_retrieved['new_name']==sample][['file_data']].values[0][0]
        length_info_file = pd_samples_retrieved[pd_samples_retrieved['new_name']==sample][['length_table']].values[0][0]
        
        print (dup_annot_file)
        print (length_info_file)
        
        ## call dup_plotter for each sample and create BioCircos plots and html entry
        
        
        
