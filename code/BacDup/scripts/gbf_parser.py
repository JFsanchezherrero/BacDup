#!/usr/bin/env python3
##############################################################
## Jose F. Sanchez & Alba Moya                              ##
## Copyright (C) 2020-2021                                  ##
##############################################################
'''
Created on 30 oct. 2020
@author: alba

Modified in March 2021
@author: Jose F. Sanchez-Herrero
'''

## useful imports
import os
import sys
import pandas as pd
import numpy as np

from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser

from HCGB.functions.aesthetics_functions import debug_message
from BacDup.scripts.functions import columns_annot_table

################################################################################     
def gbf_parser_caller(annot_file, output_path, debug):
    
    ## set output paths
    prot_file = os.path.abspath( os.path.join(output_path, 'proteins.fa'))
    csv_file = os.path.abspath( os.path.join(output_path, 'annot_df.csv'))
    csv_length = os.path.abspath( os.path.join(output_path, 'length_df.csv'))
    list_out_files = [prot_file, csv_file, csv_length]
    
    try:
        with open(prot_file, "w") as output_handle:
            SeqIO.write(
                gbf_parser(annot_file, list_out_files, debug=debug), 
                output_handle, "fasta")
    
        ## output files    
        return (list_out_files)
    except:
        return (False)
    
################################################################################
def gbf_parser(gbf_file, list_out_files, debug=False):
    
    #create an empty dataframe. 
    
    ## get common column names
    columns = columns_annot_table()
    
    annot_df = pd.DataFrame(data=None, columns=columns)
    genome_length = pd.DataFrame(data=None, columns=["length"])
    
    for rec in SeqIO.parse(gbf_file, "genbank"):
        #get genome length for BioCircos plotting  
        ID = rec.id
        genome_length.loc[ID,["length"]]=[len(rec.seq)]
        
        ## debug messages     
        if (debug):
            debug_message('GenBank record', 'yellow')
            print(rec)
            
        ## loop through features
        for feature in rec.features:
              
            #sort by CDS type. Duplicate genes analysis needs coding regions to proteins.
            if feature.type=="CDS":
                if int(feature.strand) > 0:
                    strand = "pos"
                else:
                    strand = "neg"
            
                #we create an ID for each entry     
                protID = feature.type + "_" + rec.id + "_" + str(feature.location.nofuzzy_start) + "_" + str(feature.location.nofuzzy_end) + "_" + strand
                annot_df.loc[protID, ["rec_id", "start", "end", "strand"]] = [ID, feature.location.nofuzzy_start, feature.location.nofuzzy_end, strand]
                qualif = feature.qualifiers

                ## Debug messages               
                if (debug):
                    debug_message('protID: ' + protID, 'yellow')
                    debug_message('qualif: ', 'yellow')
                    print (qualif)

                for keys, values in qualif.items():
                    #fill the dataframe info
                    if keys not in columns:
                        continue
                    annot_df.loc[protID,[keys]] = [values[0]]
                    
                    if keys=="pseudo":
                        pseudo = "True"
                        annot_df.loc[protID,["pseudo"]] = [pseudo]
            
                        #we create a FASTA file with protein sequences

                        ## Debug messages               
                        if (debug):
                            debug_message('feature: ', 'yellow')
                            print(feature)
                   
                if keys=="translation":
                    #pseudogenes have no translation item
                    gene_seq = Seq.Seq(feature.qualifiers["translation"][0])
                else:
                    pass

                yield(SeqRecord(gene_seq, protID,"",""))
                       
    ## print to file
    annot_df.to_csv(list_out_files[1], header=True)
    genome_length.to_csv(list_out_files[2], header=False)

    ## debug messages
    if (debug):
        debug_message('annot_df: ', 'yellow')
        print(annot_df)
               
    return()

############################################
def main (gbf_file, output_folder, debug=False):
    #get name
    base, ext = os.path.splitext(gbf_file)
    gbf_file = os.path.abspath(gbf_file)
    
    #create folder
    output_path = HCGB.functions.file_functions.create_folder(output_path)
    
    if (debug):
        print ("## DEBUG:")
        print ("base:" , base)
        print ("ext:" , ext)
        print ()
        
    gbf_parser_caller(gbf_file, output_path, debug)
    
############################################
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print (__doc__)
        print ("## Usage protein_gbf")
        print ("python %s gbf_file output_folder\n" %sys.argv[0])
        sys.exit()
    main(*sys.argv[1:], debug=True)