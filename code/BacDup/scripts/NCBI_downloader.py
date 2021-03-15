#usr/bin/env python
'''
This code downloads information for NCBI assembly IDs provided and updates/populates de database of interest
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''
## useful imports
import os
import sys
import pandas as pd
import ncbi_genome_download as ngd
from Bio import SeqIO
from termcolor import colored

import HCGB

###############################################################
def ngd_download(dir_path, acc_ID, data_folder, debug):
    '''
    Modified in March 2021
    @author: Jose F. Sanchez-Herrero

    Code retrieve from BacterialTyper database_generator.py script
    '''
    download = False
    print ('+ Check data for ID: ', acc_ID)
    if os.path.exists(dir_path):
        print ('+ Folder already exists: ', dir_path)
        ## get files download
        (genome, prot, gff, gbk) = get_files_download(dir_path)
        if all([genome, prot, gff, gbk]):
            download = False
        else:
            print ('+ Not all necessary data is available. Download it again.')
            download = True
    else:
        download = True
    
    if download:
        print ('+ Downloading:')
        ## download in data folder provided
        ngd.download(section='genbank', file_formats='fasta,gff,protein-fasta,genbank', assembly_accessions=acc_ID, output=data_folder, groups='bacteria')

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
        print ('+ Data is already available, no need to download it again')

###############################################################
def NCBI_download_list(strains2get, data_folder):
    ## loop through strains2get and call NCBIdownload
    for acc_ID in strains2get:
        HCGB.functions.aesthetics_functions.print_sepLine("+", 75, False)
        dir_path = data_folder + '/genbank/bacteria/' + acc_ID ## module ngd requires to download data in bacteria subfolder under genbank folder

        ## check if already exists
        if os.path.isdir(dir_path):
            print ("\n+ Data is already available in database for: " + colored(acc_ID, 'green'))

        else:
            ## download
            print ("\n+ Downloading data for: " + colored(acc_ID, 'green'))
            data_accID = NCBIdownload(acc_ID, data_folder)
            this_db = functions.get_data(data_accID, ',', 'index_col=0')
            this_db = this_db.set_index('ID')
            database_df = database_df.append(this_db)


##########################################################################################
def NCBIdownload(acc_ID, data_folder, debug):    
    '''
    Modified in March 2021
    @author: Jose F. Sanchez-Herrero

    Code retrieve from BacterialTyper database_generator.py script
    '''
    ## module ngd requires to download data in bacteria subfolder under genbank folder
    dir_path = os.path.join(data_folder, acc_ID) 
    ngd_download(dir_path, acc_ID, data_folder)
    
    ## get files download
    (genome, prot, gff, gbk) = get_files_download(dir_path)

    ## check if any plasmids downloaded
    plasmid_count = 0
    plasmid_id = []
    contig_out_file = dir_path + '/' + acc_ID + '_chromosome.fna'
    plasmid_out_file = dir_path + '/' + acc_ID + '_plasmid.fna' 
    
    ## open
    contig_out_file_handle = open(contig_out_file, 'w')
    for seq_record in SeqIO.parse(genome, "fasta"):
        plasmid_search = re.search(r".*plasmid.*", seq_record.description)
        if plasmid_search:
            ## count and get names for plasmids
            plasmid_count += 1
            name = str( seq_record.id )
            plasmid_id.append(name)
        
            ### Separate plasmids from main sequence
            plasmid_out_file_handle = open(plasmid_out_file, 'a')
            plasmid_out_file_handle.write(seq_record.format("fasta"))
            plasmid_out_file_handle.write('\n')
            plasmid_out_file_handle.close()
        else:
            contig_out_file_handle.write(seq_record.format("fasta"))
            contig_out_file_handle.write('\n')
            contig_out_file_handle.close()

    ## no plasmids found
    if plasmid_count == 0:
        plasmid_out_file = ""
        
    data2download=pd.DataFrame(columns=('ID','folder','genus','species','name','genome', 'chr', 'GFF','GBK', 'proteins','plasmids_number','plasmids_ID','plasmids'))
    data2download.loc[len(data2download)] = (acc_ID, dir_path, data.loc[acc_ID]['genus'], 
                                            data.loc[acc_ID]['species'], data.loc[acc_ID]['name'], genome, 
                                            contig_out_file, gff, prot, gbk,
                                            plasmid_count, "::".join(plasmid_id), plasmid_out_file)

    ## return data
    return(data2download)        
    
##########################################################################################
def get_files_download(folder):
    '''
    Modified in March 2021
    @author: Jose F. Sanchez-Herrero

    Code retrieve from BacterialTyper database_generator.py script
    '''
    ## check if files are gunzip
    files = os.listdir(folder)
    genome=""
    prot=""
    gff=""
    gbk=""
    for f in files:
        if f.endswith('genomic.fna'):
            genome = os.path.join(folder, f)
        elif f.endswith('genomic.gff'):
            gff = os.path.join(folder, f)
        elif f.endswith('genomic.gbk'):
            gbk = os.path.join(folder, f)
        elif f.endswith('genomic.gbff'):
            gbk = os.path.join(folder, f)
        elif f.endswith('protein.faa'):
            prot = os.path.join(folder, f)

    return(genome, prot, gff, gbk)           


###############################################################
def help_options():
    print ("\nUSAGE: python %s ID_file folder\n"  %os.path.realpath(__file__))
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
    NCBI_download_list(strains2get, folder)
    
    print ("+ Data has been retrieved: \n", data)


############################################################
if __name__== "__main__":
    main()
    
    
## Download all genomes from a taxa and descendent
## https://www.biostars.org/p/302533/    
