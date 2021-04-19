#!/usr/bin/env python3
##############################################################
## Jose F. Sanchez & Alba Moya                              ##
## Copyright (C) 2020-2021                                  ##
##############################################################
'''
Created on 28 oct. 2020
@author: alba

Modified in March 2021
@author: Jose F. Sanchez-Herrero
'''

## useful imports
import sys
import os
import pandas as pd
import numpy as np
import HCGB

from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
from BCBio import GFF
from BacDup.scripts.functions import columns_annot_table

   
##################################################
def gff_parser_caller(gff_file, ref_file, output_path, debug):
    '''This function calls the actual gff parser
    
    It serves as the entry point either from a module or system call
    '''
    
    ## set output paths
    prot_file = os.path.abspath( os.path.join(output_path, 'proteins.fa'))
    csv_file = os.path.abspath( os.path.join(output_path, 'annot_df.csv'))
    csv_length = os.path.abspath( os.path.join(output_path, 'length_df.csv'))
    list_out_files = [prot_file, csv_file, csv_length]
    
    with open (ref_file) as in_handle:
        ref_recs = SeqIO.to_dict(SeqIO.parse(in_handle, "fasta"))
    
    ## debug messages     
    if (debug):
        debug_message('GenBank record', 'yellow')
        print (ref_recs)

    ## parse
    with open(prot_file, "w") as out_handle:
        SeqIO.write(protein_recs(gff_file, ref_recs, list_out_files, debug=debug), out_handle, "fasta")
    
    ## return information
    return (list_out_files)

############################################################       
def protein_recs(gff_file, ref_recs, list_out_files, debug=False):
    '''GFF parser to retrieve proteins and annotation
    '''
    
    #create an empty dataframe. 
    columns = columns_annot_table()    ## get common column names
    annot_df = pd.DataFrame(data=None, columns=columns)
    genome_length = pd.DataFrame(data=None, columns=["length"])
    
    with open(gff_file) as in_handle:
        ##parse the output. Generate SeqRecord and SeqFeatures for predictions
        ##sort by CDS type. Duplicate genes analysis just needs coding regions to proteins.
        limit_info = dict(gff_type=["CDS"])    
        for rec in GFF.parse(in_handle, limit_info = limit_info, base_dict=ref_recs):
            #get genome length for BioCircos plotting  
            ID = rec.id
            genome_length.loc[ID,["length"]]=[len(rec.seq)]
    
            ## debug messages     
            if (debug):
                debug_message('GFF record', 'yellow')
                print(rec)
           
            for feature in rec.features:
                ## Debug messages               
                if (debug):
                    debug_message('feature: ', 'yellow')
                    print(feature)

                ## strand
                if feature.strand == -1:
                    strand = "neg"
                else:   
                    strand = "pos"
                    
                #we create an ID for each entry     
                protID = feature.type + "_" + rec.id + "_" + str(feature.location.nofuzzy_start) + "_" + str(feature.location.nofuzzy_end) + "_" + strand
                annot_df.loc[protID, ["rec_id", "type", "start", "end", "strand"]] = [ID, feature.type, feature.location.nofuzzy_start, feature.location.nofuzzy_end, strand]
                qualif = feature.qualifiers
                ## Debug messages               
                if (debug):
                    debug_message('protID: ' + protID, 'yellow')
                    debug_message('qualif: ', 'yellow')
                    print (qualif)

                ## loop
                for keys, values in qualif.items():
                    #fill the dataframe info
                    if keys == "Note":
                        continue
                    annot_df.loc[protID,[keys]] = [values[0]]
                    
                ## get gene sequence
                gene_seq = Seq.Seq(str(rec.seq[feature.location.nofuzzy_start:feature.location.nofuzzy_end]))

                ## Debug messages               
                if (debug):
                    debug_message('gene_seq: ' + protID, 'yellow')
                    print (gene_seq)

                if feature.type == "CDS":
                    if feature.strand == -1:
                        gene_seq = gene_seq.reverse_complement()
                    
                    #delete STOP symbols
                    protein_seq = gene_seq.translate(to_stop=True)
                    
                    yield(SeqRecord(protein_seq, protID, "", ""))

    ## print additional information
    annot_df.to_csv(list_out_files[1], header=True)
    genome_length.to_csv(list_out_files[2], header=True)
    
    #get genome length for BioCircos plotting  
    #genome_length = pd.DataFrame(data=None, columns=["length"])
    #ID = rec.id
    #length = len(rec.seq)
    #genome_length.loc[ID,["length"]]=[length]
    #csv_length = "%s/%s_length.csv" % (output_path, rec.id)            
    #genome_length.to_csv(csv_length, header=True)
    
    ## debug messages
    if (debug):
        debug_message('annot_df: ', 'yellow')
        print(annot_df)
    
    ## empty return
    return()


#################################################################
def main (gff_file, ref_file, output_folder, debug=False):
    #get name
    base, ext = os.path.splitext(gff_file)
    gff_file = os.path.abspath(gff_file)
    
    #create folder
    output_path = HCGB.functions.file_functions.create_folder(output_path)
    
    if (debug):
        print ("## DEBUG:")
        print ("base:" , base)
        print ("ext:" , ext)
        print ()
        
    gff_parser_caller(gff_file, ref_file, output_path, debug)
    

################################################################################
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print (__doc__)
        
        print ("## Usage gff_parser")
        print ("python %s gff_file ref_fasta_file output_folder\n" %sys.argv[0])

        sys.exit()

    main(*sys.argv[1:], debug=True)
    #main(*sys.argv[1:])

    # la variable debug no es obligatoria. tiene un "por defecto definido"
    # Se utiliza el "=" para indicar el default.
        
        
        
