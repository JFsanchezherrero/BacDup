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
## useful imports
import os
import sys
import argparse
import time
from Bio import SeqIO
import HCGB
from HCGB.functions.aesthetics_functions import debug_message
import HCGB.functions.time_functions as time_functions
from termcolor import colored
import pandas as pd

## my modules
import BacDup
import BacDup.scripts.gbf_parser as gbf_parser
import BacDup.scripts.gff_parser as gff_parser
import BacDup.scripts.format_checker as format_checker
import BacDup.scripts.functions as BacDup_functions
import BacDup.scripts.taxonomy_retrieval as taxonomy_retrieval

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
    time_functions.print_time()
    
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

    ## time stamp
    start_time_partial = time_functions.timestamp(start_time_total)
    
    #################################################
    ## Parse and obtain the type of input information provided
    #################################################
    df_accID = parse_options(arg_dict)
    ## pd.DataFrame: 'new_name','folder','genus',
    ##               'species','taxonomy','genome', 
    ##               'annot_file','format_annot_file', 'proteins',
    ##               'plasmids_number','plasmids_ID'))
    
    ## time stamp
    start_time_partial = time_functions.timestamp(start_time_partial)
    
    ## parse information accordingly
    parse_information(arg_dict, df_accID, outdir)

    ### report generation
    HCGB.functions.aesthetics_functions.boxymcboxface("Summarizing input files")
    outdir_report = HCGB.functions.files_functions.create_subfolder("report", outdir)

    input_report = HCGB.functions.files_functions.create_subfolder("input", outdir_report)
    
    ## add df_accID.loc[sample,] information as csv into input folder
    df_accID.to_csv(os.path.join(input_report, 'info.csv'), index=True, header=True)
    
    ## maybe add a summary of the files?
    
    print ("\n*************** Finish *******************")
    start_time_partial = time_functions.timestamp(start_time_total)

    print ("+ Exiting Input module.")
    return()

################################################################################
def parse_information(arg_dict, df_accID, outdir):

    ### Parse df_accID
    dict_input_folders = HCGB.functions.files_functions.outdir_project(outdir, arg_dict.project, df_accID, "input", arg_dict.debug)
    dict_parse_folders = HCGB.functions.files_functions.outdir_project(outdir, arg_dict.project, df_accID, "parse", arg_dict.debug)

    ## debug messages
    if (arg_dict.debug):
        debug_message('+++++++++++++++++++++++++++++++')
        print("dict_input_folders")
        print(dict_input_folders)
        print("dict_parse_folders")
        print(dict_parse_folders)

    ## parse each sample retrieved
    for sample, folder_input in dict_input_folders.items():

        if (arg_dict.debug):
            debug_message('sample: ' + sample, 'yellow')
            debug_message('folder_input: ' + folder_input, 'yellow')
            debug_message('folder_parse: ' + dict_parse_folders[sample], 'yellow')
            debug_message('annot_file: ' + df_accID.loc[sample, 'annot_file'], 'yellow')
            debug_message('genome' + df_accID.loc[sample, 'genome'], 'yellow')

        ## timestamps 
        input_timestamp = os.path.join(folder_input, '.success')
        parse_timestamp = os.path.join(dict_parse_folders[sample], '.success')
        
        print()
        print ("\t+ Parsing sample: " + sample)
        
        if (not HCGB.functions.files_functions.is_non_zero_file(parse_timestamp) and not HCGB.functions.files_functions.is_non_zero_file(input_timestamp)):
        
            ## TODO: Set threads to use in parallel
            process_OK = parse_annot_file(sample, folder_input, df_accID.loc[sample, 'annot_file'], dict_parse_folders[sample], arg_dict.debug, df_accID.loc[sample, 'genome'])
            
            if (process_OK):
            
                ## link or copy annotation file into folder_input
                HCGB.functions.files_functions.get_symbolic_link_file(df_accID.loc[sample, 'annot_file'], folder_input)
                
                ## add df_accID.loc[sample,] information as csv into input folder
                df_accID.loc[sample,].to_csv(os.path.join(folder_input, 'info.csv'), index=True, header=True)
                
                ## print time stamp
                time_functions.print_time_stamp(input_timestamp)
        
                ## print time stamp
                time_functions.print_time_stamp(parse_timestamp)
            else:
                print(colored("\t+ Some error occurred for sample %s while parsing input options" %sample, 'red'))
                
                ## print time stamp
                time_functions.print_time_stamp(os.path.join(folder_input, '.fail'))
        
                ## print time stamp
                time_functions.print_time_stamp(os.path.join(dict_parse_folders[sample], '.fail'))
        else:
            read_time = time_functions.read_time_stamp(parse_timestamp)
            print (colored("\t+ Input parsing already available for sample %s [%s]" %(sample, read_time), 'green'))
            print()

##########################
def input_help():
    # TODO: Update this information
    print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))
    
    
    compt = {}
    compt["fasta"] = [".fa", ".faa", ".mpfa", ".fna", ".fsa", ".fas", ".fasta"]
    compt["genbank"] = [".genbank", ".gb", ".gbf", ".gbff", ".gbk"]
    compt["GFF"] = [".gff"]
    
    
    
#############################
def parse_annot_file(name, folder_out_input, annot_file, output_path, Debug, ref_file=""):
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
    
    ## check file integrity: exists & non-zero
    if (BacDup_functions.file_readable_check(annot_file)):
        ## check format; call parser
        format = format_checker.is_format(annot_file, Debug)
        
        ## debug messages
        if (Debug):
            debug_message('\nformat_checker.is_format function call:', 'yellow')
            debug_message('format: ' + format, 'yellow')
            
        ## parse gbk or gff        
        if (format=='gbk'):
            print (colored('\t* GenBank format file:........[OK]', 'green'))
            
            ## TODO: print details available within GenBank:
            # Accession, Bioproject,
            # Reference, Authors, Title, Journal,
            # Comment
            
            return(gbf_parser.gbf_parser_caller(annot_file, output_path, Debug))
        
        elif(format=='gff'):
            print (colored('\t* GFF format file:.......[OK]', 'green'))
            if (HCGB.functions.files_functions.is_non_zero_file(ref_file)):
                return(gff_parser.gff_parser_caller(annot_file, ref_file, output_path, Debug))
            else:
                print(colored("ERROR: No genome reference file provided for this GFF annotation. Check input options provided.","red"))
                exit()
        
        ## not valid via this option
        else:
            print(colored("ERROR: not valid via this option","red"))
            exit()
            
    ## not accessible for this sample
    else:
        return (False)
    
####################################
def parse_options(arg_dict):
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
            
        ## create dataframe df_accID to match other formats
        df_accID=pd.DataFrame(columns=(BacDup_functions.columns_accID_table()))
        
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
            taxonomy = "" 
            organism  = ""
            taxonomy_string = ""
            genus = ""
            if (format == 'gbk'):
                ## get information from each sample
                (taxonomy, organism) = BacDup.scripts.functions.get_gbk_information(file_annot, arg_dict.debug)
                ## plasmid_count, plasmid_id not available
                
            elif (format == 'gff'):
                if (arg_dict.ref_file):
                    arg_dict.ref_file = os.path.abspath(arg_dict.ref_file)
                    BacDup_functions.file_readable_check(arg_dict.ref_file)

                    if (arg_dict.batch):
                        ref_entries = HCGB.functions.main_functions.file2dictionary(arg_dict.ref_file, ',')
                        genome = ref_entries[name]
                    else:
                        genome = arg_dict.ref_file

            ## save into dataframe
            if len(taxonomy) > 1:
                genus = taxonomy[-1]
                taxonomy_string = ";".join(taxonomy)
                
            dir_path = os.path.abspath(os.path.dirname(file_annot))
            df_accID.loc[len(df_accID)] = (name, dir_path, genus, organism, taxonomy_string, genome, 
                                           file_annot, format, prot, 
                                           plasmid_count, ";".join(plasmid_id))

    ## --------------------------------------- ##
    ## NCBI RefSeq/Genbank IDs: GCA_XXXXXXXX.1; GCF_XXXXXXXXX.1
    ## --------------------------------------- ##
    elif (arg_dict.GenBank_id):
        ## get database path
        if (arg_dict.db_folder):
            db_folder = HCGB.functions.files_functions.create_folder(os.path.abspath(arg_dict.db_folder))
        else:
            db_folder = HCGB.functions.files_functions.create_subfolder("db", os.path.abspath(arg_dict.output_folder))

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
            strains2get = list(filter(None, strains2get))
                
            ## debug messages
            if (arg_dict.debug):
                debug_message('strains2get: ' + str(strains2get), 'yellow')
                
            ## call NCBI_downloader
            df_accID = BacDup.scripts.NCBI_downloader.NCBI_download_list(strains2get, db_folder, arg_dict.debug, arg_dict.assembly_level)
            
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
        
        ## assume all belong to same superkingdom if children of same tax_id
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
        (strains2get, allstrains_available) = taxonomy_retrieval.get_GenBank_ids(db_folder, string_info_total, int(arg_dict.k_random), 
                                                          arg_dict.debug, assembly_level_given=arg_dict.assembly_level,
                                                          group_given=group_obtained, section_given=arg_dict.section)

        ## print list and dictionary of possible and selected taxIDs
        outdir = os.path.abspath(arg_dict.output_folder)
        final_dir = HCGB.functions.files_functions.create_subfolder("info", outdir)
        input_info_dir = HCGB.functions.files_functions.create_subfolder("input", outdir)
        HCGB.functions.main_functions.printList2file(os.path.join(input_info_dir, 'Downloaded.txt'), strains2get)
        HCGB.functions.main_functions.printList2file(os.path.join(input_info_dir, 'all_entries.txt'), allstrains_available)
        
        ## save into file
        file_info = os.path.join(input_info_dir, 'info.txt')
        
        #################
        ## call NCBI_downloader
        #################
        df_accID = BacDup.scripts.NCBI_downloader.NCBI_download_list(strains2get, db_folder, arg_dict.debug, arg_dict.assembly_level)

    ## --------------------------------------- ##
    ## Previous BacDup analysis folder
    ## --------------------------------------- ##
    ## TODO
    elif (arg_dict.project):
        print (colored('\t* A previous BacDup analysis project folder:.......[OK]', 'green'))
        ## create df_accID to store data
        ## TODO

    ## Returns dataframe with information
    
    df_accID = df_accID.set_index('new_name')
    return (df_accID)