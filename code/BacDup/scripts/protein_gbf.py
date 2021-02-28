'''
Created on 30 oct. 2020

@author: alba
'''

import os
import sys
import pandas as pd
import numpy as np
import output

import Bio
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser

#####
def main (gbf_file, output_folder, debug=False):
    #get name
    base, ext = os.path.splitext(gbf_file)
    gbf_file = os.path.abspath(gbf_file)
    
    #create folder
    output_path = os.path.abspath(output_folder)
    output.create_folder(output_path)
    
    if (debug):
        print ("## DEBUG:")
        print ("base:" , base)
        print ("ext:" , ext)
        print ()
        
    gbf_parser_caller(gbf_file, output_path, debug)
    
 
 
#####     
def gbf_parser_caller(gbf_file, output_path, debug):
    
    prot_file = "%s/proteins.fa" % output_path
    csv_file = "%s/df.csv" % output_path            

    
    with open(prot_file, "w") as output_handle:
        SeqIO.write(gbf_parser(gbf_file, output_path, debug=debug), output_handle, "fasta")
   
        
    return (prot_file, csv_file)
    
#####
def gbf_parser(gbf_file, output_path, debug=False):
    
    #create an empty dataframe. The most important info as first columns
    #all the other file info will be filled next
    columns = ['rec_id', 'locus_tag', 'protein_id', 'gene', 'start',
               'end', 'strand', 'pseudo', 'product', 'db_xref',
               'EC_number', 'old_locus_tag', 'inference']
    annot_df = pd.DataFrame(data=None, columns=columns)
    genome_length = pd.DataFrame(data=None, columns=["length"])
    
    for rec in SeqIO.parse(gbf_file, "genbank"):
        #get genome length for BioCircos plotting  
        ID = rec.id
        genome_length.loc[ID,["length"]]=[len(rec.seq)]
                
        if (debug):
            print("## DEBUG: rec")
            print(rec)
            print()
            
        for feature in rec.features:
              
            #sort by CDS type. Duplicate genes analysis needs coding regions to proteins.
            if feature.type=="CDS":
                if int(feature.strand) > 0:
                    strand = "pos"
                else:
                    strand = "neg"
            
                #we create an ID for each entry     
                protID = feature.type + "_" + rec.id + "_" + str(feature.location.nofuzzy_start) + "_" + str(feature.location.nofuzzy_end) + "_" + strand
#                 if (debug):
#                     print("#DEBUG: protID")
#                     print(protID)
#                     print()
                        
                annot_df.loc[protID, ["rec_id", "start", "end", "strand"]] = [ID, feature.location.nofuzzy_start, feature.location.nofuzzy_end, strand]
                qualif = feature.qualifiers
#                 if (debug):
#                     print("#DEBUG: Qualifiers")
#                     print(qualif)
#                     print()
                for keys, values in qualif.items():
                 
                        #fill the dataframe info
                    if keys not in columns:
                        continue
                    annot_df.loc[protID,[keys]] = [values[0]]
                    
                    if keys=="pseudo":
                        pseudo = "True"
                        annot_df.loc[protID,["pseudo"]] = [pseudo]
            
                        #we create a FASTA file with protein sequences
            
#                     if (debug):
#                         print ("## DEBUG: feature")
#                         print (feature)
#                         print ()
                    
#                         print (feature.location)
#                         print (feature.location.nofuzzy_start)
#                         print (feature.location.nofuzzy_end)
#                         print (feature.strand)
#                         print (feature.id)
#                         print (feature.qualifiers)   
                   
                if keys=="translation":
                    #pseudogenes have no translation item
                    gene_seq = Seq.Seq(feature.qualifiers["translation"][0])
                else:
                    pass
                yield(SeqRecord(gene_seq, protID,"",""))
                       
    
    csv_file = "%s/df.csv" % output_path            
    annot_df.to_csv(csv_file, header=True)
    csv_length = "%s/length.csv" % (output_path)            
    genome_length.to_csv(csv_length, header=False)
    
    
  
    if (debug):
        print("## DEBUG: dataframe")
        print(annot_df)
               
      
  
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print (__doc__)
        
       
        print ("## Usage protein_gbf")
        print ("python %s gbf_file output_folder\n" %sys.argv[0])
       

        sys.exit()

    main(*sys.argv[1:], debug=True)
  
    
    
    

            
            
