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
from Bio.SeqFeature import FeatureLocation, ExactPosition
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
    
    ## create dataframe. 
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
                genome_seq = rec.seq[feature.location.nofuzzy_start:feature.location.nofuzzy_end]
                if int(feature.strand) > 0:
                    strand = "pos"
                else:
                    strand = "neg"
                    genome_seq = genome_seq.reverse_complement()
                    
                #we create an ID for each entry     
                protID = feature.type + "_" + rec.id + "_" + str(feature.location.nofuzzy_start) + "_" + str(feature.location.nofuzzy_end) + "_" + strand
                annot_df.loc[protID, ["rec_id", "start", 
                                      "end", "strand"]] = [ID, feature.location.nofuzzy_start, 
                                                           feature.location.nofuzzy_end, strand]
                qualif = feature.qualifiers
                pseudo=False
                
                ## Debug messages               
                if (debug):
                    debug_message('protID: ' + protID, 'yellow')
                    
                    debug_message('qualif: ', 'yellow')
                    print (qualif)
                    
                    debug_message('feature: ', 'yellow')
                    print(feature)
                    
                    debug_message('genome_seq: ', 'yellow')
                    print(genome_seq)
                
                pseudo_seq = ""
                ## fill datafarme
                for keys, values in qualif.items():
                    if keys not in columns:
                        continue
                    
                    ## Save keys into dataframe
                    annot_df.loc[protID,[keys]] = [values[0]]
                    
                    ####################################
                    ## Pseudogenes:
                    ####################################
                    if keys=="pseudo":
                        pseudo = True
                        
                        ## set pseudo True/False
                        annot_df.loc[protID,["pseudo"]] = ["True"]            
                        table_code = feature.qualifiers["transl_table"][0]
                        pseudo_seq = genome_seq.translate(table=table_code, to_stop=False)
                        
                        ## Debug messages               
                        if (debug):
                            print("***************************************")
                            debug_message('Pseudogene: ', 'yellow')
                            print("***************************************")
                            debug_message('feature.location.nofuzzy_start: ', 'yellow')
                            print(feature.location.nofuzzy_start)
                            debug_message('feature.location.nofuzzy_end: ', 'yellow')
                            print(feature.location.nofuzzy_end)
                            debug_message('Translation table code: ', 'yellow')
                            print(table_code)
                            debug_message('genome_seq: ', 'yellow')
                            print(genome_seq)
                            debug_message('pseudo_seq: ', 'yellow')
                            print(pseudo_seq)
                            
                   
                ## create a sequence fasta entry
                if (pseudo):
                    # Pseudogenes have no translation item
                    # set translated CDS even including *
                    if len(pseudo_seq)!=0:
                        gene_seq = pseudo_seq
                    else:
                        ## sometimes it might fail
                        gene_seq=Seq.Seq('***')

                else:
                    ## CDS provided by genbank
                    gene_seq = Seq.Seq(feature.qualifiers["translation"][0]) 
                
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