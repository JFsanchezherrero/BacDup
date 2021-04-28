#!/usr/bin/env python3
'''
Created and modified in March 2021
@author: Jose F. Sanchez-Herrero

Original code retrieved from BacterialTyper database_generator.py script.

This code downloads information for NCBI assembly IDs provided
'''
## useful imports
import os
import re
import sys
import pandas as pd
import ncbi_genome_download as ngd
from termcolor import colored

import BacDup
import HCGB
from HCGB.functions.aesthetics_functions import debug_message

###############################################################
def get_section(acc_ID, debug):
    if (acc_ID.startswith("GCA")):
        return('genbank')
    else:
        return('refseq')

###############################################################
def ngd_download(section_given, acc_ID, data_folder, debug, section='genbank', assembly_level='complete', group_given='bacteria'):
    '''
    Function that calls and retrieves data from NCBI using python package ngd.
    
    :param acc_ID:
    :param data_folder: Folder to store data. 
    :param debug: True/false for debugging messages
    
    :attention Module ngd requires to download data in bacteria/archaea subfolder under genbank or refseq folder.
    '''
    ##################################
    ## check if necessary to download
    ##################################
    
    ## get path
    print ('+ Check data for ID: ', acc_ID)
    dir_path = os.path.join(data_folder, section_given, group_given, acc_ID)
    
    ## check if previously download
    download = False
    if os.path.exists(dir_path):
        print ('+ Folder already exists: ', dir_path)
        ## get files download
        (genome, prot, gff, gbk) = BacDup.scripts.functions.get_files_annotation(dir_path, debug)
        if (gbk): ## Only genbank format file is required
            download = False
        else:
            print ('+ Not all necessary data is available. Download it again.')
            download = True
    else:
        download = True
    
    ## download data
    if download:
        print ("\n+ Downloading data for: " + colored(acc_ID, 'green'))

        ## download in data folder provided
        if (debug):
            debug_message("ngd.download call", color="yellow")
            debug_message("dir_path: " + dir_path, color="yellow")
            debug_message("section_given: " + section_given, color="yellow")
            
        ## download
        if debug:
            debug_message("section='%s', file_formats='genbank', assembly_level=%s, assembly_accessions=%s, output=%s, groups=%s" %(section_given, assembly_level, acc_ID, data_folder, group_given), color="yellow")
        
        try:
            ngd.download(section=section_given, file_formats='genbank', 
                     assembly_levels=assembly_level,
                     assembly_accessions=acc_ID, output=data_folder, groups=group_given)
        except:
            raise ("A problem occurred when contacting NCBI for downloading id (%s) from %s" %(acc_ID, section_given))
            
        
        ## return empty
        if not os.path.isdir(dir_path):
            return False
    
        ## check if files are gunzip
        files = os.listdir(dir_path)
        files_list = []        
        for f in files:
            if f.endswith('gz'):
                files_list.append(f)
                print ("\t- Extracting files: ", f)
                HCGB.functions.files_functions.extract(dir_path + '/' + f, dir_path)
                #os.remove(dir_path + '/' + f)
    
    ## skip
    else:
        print ('\t+ Data is already available, no need to download it again')
    
    print()
    ## return path where data is
    return (dir_path)

###############################################################
def NCBI_download_list(strains2get, data_folder, Debug, assembly_level):
    '''
    Function to call ngd_download function given a list of IDs.
    Returns dataframe containing all samples downloaded.
    '''
    
    ## debug messages
    if Debug:
        debug_message("******************************************")
        debug_message("NCBI_download_list function call", color="yellow")
        debug_message("strains2get", color="yellow")
        debug_message(strains2get, color="yellow")
        debug_message("data_folder", color="yellow")
        debug_message(data_folder, color="yellow")
        debug_message("******************************************")
    
    ## get unique values
    strains2get = list(set(strains2get))
    
    ## loop through strains2get and call NCBIdownload
    database_df = pd.DataFrame()
    for acc_ID in strains2get:
        ## TODO: Set threads to use in parallel
        HCGB.functions.aesthetics_functions.print_sepLine("+", 75, False)
        data_accID = NCBIdownload(acc_ID, data_folder, Debug, assembly_level)
        database_df = database_df.append(data_accID)

    ## debug messages
    if Debug:
        debug_message("******************************************")
        debug_message("Return info: NCBI_download_list", color="yellow")
        debug_message("database_df", color="yellow")
        HCGB.functions.main_functions.print_all_pandaDF(database_df)
        debug_message("******************************************")
        
    ## return dataframe containing all information
    return(database_df)

##########################################################################################
def NCBIdownload(acc_ID, data_folder, debug, assembly_level='complete', group='bacteria'):    
    '''
    Code retrieve from BacterialTyper database_generator.py script
    '''
    ## module ngd requires to download data in bacteria/archaea subfolder under genbank/refseq folder
    ## get details
    section_given = get_section(acc_ID, debug)
    group_accID = NCBI_get_info_GenbankID(data_folder, acc_ID, debug)
    
    if debug:
        debug_message("-----------------------------------------")
        debug_message("NCBIdownload function call", color="yellow")
        debug_message("section_given: " + section_given, color="yellow")
        debug_message("group_accID: " + group_accID, color="yellow")

    ## save into dataframe
    dataDownloaded=pd.DataFrame(columns=(BacDup.scripts.functions.columns_accID_table()))
    
    ## download accID
    dir_path = ngd_download(section_given, acc_ID, data_folder, debug, assembly_level=assembly_level, group_given=group_accID)
    
    ## return empty
    if not os.path.isdir(dir_path):
        dataDownloaded.loc[len(dataDownloaded)] = (acc_ID, '', '', '', '', '', '', '', '', '', '')
        return(dataDownloaded)
    
    ## get files download
    (genome, prot, gff, gbk) = BacDup.scripts.functions.get_files_annotation(dir_path, debug)

    ## check if any plasmids downloaded
    (plasmid_count, plasmid_id) = BacDup.scripts.functions.get_plasmids(genome, debug)
        
    ## get information from each sample
    (taxonomy, organism) = BacDup.scripts.functions.get_gbk_information(gbk, debug)
    taxonomy_string = ";".join(taxonomy)

    ## save into dataframe
    dataDownloaded.loc[len(dataDownloaded)] = (acc_ID, dir_path, taxonomy[-1], organism, taxonomy_string, genome, gbk, 'gbk', prot, plasmid_count, ";".join(plasmid_id))

    ## return dataframe containing all information
    if debug:
        debug_message("-----------------------------------------")
        debug_message("Return info NCBIdownload", color="yellow")
        debug_message("dataDownloaded", color="yellow")
        debug_message(dataDownloaded, color="yellow")
        debug_message("-----------------------------------------")
        
    return(dataDownloaded)        
    
###############################################################
def NCBI_get_info(section_given, data_folder, tax_ID_list, debug, assembly_level_given ='complete', group_given='bacteria'):
    '''This function uses ncbi_genome_download to 
    create a dry run and return information of each entry provided
    '''
    ## import module and class
    import ncbi_genome_download
    from ncbi_genome_download.config import NgdConfig
    
    try:
        ngd_config = NgdConfig.from_kwargs(section=section_given, 
                     file_formats='genbank',
                     taxids=tax_ID_list,
                     output=data_folder,
                     dry_run=True, 
                     assembly_levels=assembly_level_given,
                     groups=group_given)
        info = ncbi_genome_download.core.select_candidates(ngd_config)
        
    except:
        raise "**** ERROR: Something happen while connecting to NCBI... ***"
        exit()
        return (False)

    ####
    if (len(info)) < 1:
        print (colored("No entries matched your filter. Please check the input options provided", 'yellow'))
        exit()

    ## fill dictionary to simplify
    dict_entries = {}
    for entry, _ in info:
        strain_name = ncbi_genome_download.core.get_strain(entry)
        ## debug messagess
        if debug:
            debug_message("", 'yellow')
            print(entry)
            string = entry['assembly_accession'] + '\t' + entry['organism_name'] + '\t' + strain_name
            debug_message(string, 'yellow')
            debug_message(".....................................................................\n", 'yellow')
            
        ## fill dictionary
        dict_entries[entry['assembly_accession']] = (entry['organism_name'], strain_name)
        
    ## return
    return(dict_entries)

###############################################################
def NCBI_get_info_GenbankID(data_folder, acc_ID, debug, assembly_level_given ='complete'):

    section_given = get_section(acc_ID, debug)

    if debug:
        debug_message("-----------------------------------------")
        debug_message("NCBI_get_info_GenbankID function call", color="yellow")


    ## import module and class
    import ncbi_genome_download
    from ncbi_genome_download.config import NgdConfig
    tries = ['bacteria', 'archaea']
    for entry_tried in tries:
        if debug:
            debug_message("Trying with: " + entry_tried, color="yellow")

        ngd_config = NgdConfig.from_kwargs(section=section_given, 
                     file_formats='genbank',
                     assembly_accessions=acc_ID,
                     output=data_folder,
                     dry_run=True, 
                     groups=entry_tried)
        info = ncbi_genome_download.core.select_candidates(ngd_config)
        if info:
            debug_message("It worked!", color="yellow")
            return(entry_tried)
        
    raise "**** ERROR: Something happen while connecting to NCBI... ***"
    exit()
    return (False)

###############################################################
def help_options():
    print ("\nUSAGE: python %s ID_file folder Debug[True/False]\n"  %os.path.realpath(__file__))
    print ("---------------")
    print ("ID_file example: list of Genbank or Refseq IDs")    
    print ("---------------")
    print ("GCx_0000001")
    print ("GCx_0000002")
    print ("GCx_0000003")
    print ("GCx_0000004")
    print ("GCx_0000005")
    print ("*******")
    print ("______________")
        
############################################################
def main():
    ## this code runs when call as a single script

    ## control if options provided or help
    if len(sys.argv) > 1:
        print ("")
    else:
        help_options()
        exit()
    
    ## get arguments provided
    ID_file = os.path.abspath(sys.argv[1])
    folder = os.path.abspath(sys.argv[2])
    
    ## get file information
    strains2get = HCGB.functions.main_functions.readList_fromFile(ID_file)
    
    ## debug messages
    Debug=False
    if (sys.argv[3]=="True"):
        print ('*******************************')
        debug_message("Mode ON")
        print ('*******************************')
        Debug=True
    
    data = NCBI_download_list(strains2get, folder, Debug)
    print ("+ Data has been retrieved.\n")
    

############################################################
if __name__== "__main__":
    main()
    
    
## Download all genomes from a taxa and descendent
## https://www.biostars.org/p/302533/    
