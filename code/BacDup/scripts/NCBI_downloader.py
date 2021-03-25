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
def ngd_download(acc_ID, data_folder, debug):
    '''
    Function that calls and retrieves data from NCBI using python package ngd.
    
    :param acc_ID:
    :param data_folder: Folder to store data. 
    :param debug: True/false for debugging messages
    
    :attention Module ngd requires to download data in bacteria subfolder under genbank or refseq folder.
    '''

    ## check if necessary to download
    download = False
    print ('+ Check data for ID: ', acc_ID)
    if (acc_ID.startswith("GCA")):
        dir_path = os.path.join(data_folder, 'genbank', 'bacteria', acc_ID)
    else:
        dir_path = os.path.join(data_folder, 'refseq', 'bacteria', acc_ID)
    
    
    if os.path.exists(dir_path):
        print ('+ Folder already exists: ', dir_path)
        ## get files download
        (genome, prot, gff, gbk) = BacDup.scripts.functions.get_files_annotation(dir_path, debug)
        if all([genome, prot, gff, gbk]):
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
            
        ## download
        
        ## check if genbank (GCA_.*) or refseq ID (GCF_.*)
        if (acc_ID.startswith("GCA")):
            if debug:
                debug_message("section='genbank', file_formats='fasta,gff,protein-fasta,genbank', assembly_accessions=%s, output=%s, groups='bacteria'" %(acc_ID, data_folder), color="yellow")
            
            ngd.download(section='genbank', file_formats='fasta,gff,protein-fasta,genbank', 
                         assembly_accessions=acc_ID, output=data_folder, groups='bacteria')
        else:
            if debug:
                debug_message("section='refseq', file_formats='fasta,gff,protein-fasta,genbank', assembly_accessions=%s, output=%s, groups='bacteria'" %(acc_ID, data_folder), color="yellow")
        
            ngd.download(section='refseq', file_formats='fasta,gff,protein-fasta,genbank', 
                         assembly_accessions=acc_ID, output=data_folder, groups='bacteria')
        
        ## check if files are gunzip
        files = os.listdir(dir_path)
        files_list = []        
        for f in files:
            if f.endswith('gz'):
                files_list.append(f)
                print ("\t- Extracting files: ", f)
                HCGB.functions.files_functions.extract(dir_path + '/' + f, dir_path)
                #os.remove(dir_path + '/' + f)
    else:
        print ('\t+ Data is already available, no need to download it again')
    
    print()
    ## return path where data is
    return (dir_path)

###############################################################
def NCBI_download_list(strains2get, data_folder, Debug):
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
        data_accID = NCBIdownload(acc_ID, data_folder, Debug)
        data_accID = data_accID.set_index('new_name')
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
def NCBIdownload(acc_ID, data_folder, debug):    
    '''
    Code retrieve from BacterialTyper database_generator.py script
    '''
    ## module ngd requires to download data in bacteria subfolder under genbank folder
    dir_path = ngd_download(acc_ID, data_folder, debug)
    
    ## get files download
    (genome, prot, gff, gbk) = BacDup.scripts.functions.get_files_annotation(dir_path, debug)

    ## check if any plasmids downloaded
    (plasmid_count, plasmid_id) = BacDup.scripts.functions.get_plasmids(genome, debug)
        
    ## get information from each sample
    (taxonomy, organism) = BacDup.scripts.functions.get_gbk_information(gbk, debug)
    taxonomy_string = ";".join(taxonomy)

    ## save into dataframe
    dataDownloaded=pd.DataFrame(columns=('new_name','folder','genus','species','taxonomy','genome', 'GFF','GBK', 'proteins','plasmids_number','plasmids_ID'))
    dataDownloaded.loc[len(dataDownloaded)] = (acc_ID, dir_path, taxonomy[-1], organism, taxonomy_string, genome, gff, prot, gbk, plasmid_count, ";".join(plasmid_id))

    ## return dataframe containing all information
    if debug:
        debug_message("-----------------------------------------")
        debug_message("Return info NCBIdownload", color="yellow")
        debug_message("dataDownloaded", color="yellow")
        debug_message(dataDownloaded, color="yellow")
        debug_message("-----------------------------------------")
        
    return(dataDownloaded)        
    

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
