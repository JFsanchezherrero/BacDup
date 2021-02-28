'''
Created on 28 oct. 2020

@author: alba
'''
import sys
import os
import pandas as pd
import numpy as np
import output
import Bio
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
import BCBio
from BCBio import GFF

#### take care of the inferred parents ####
#####
def main (gff_file, ref_file, output_folder, debug=False):
    #get name
    base, ext = os.path.splitext(gff_file)
    gff_file = os.path.abspath(gff_file)
    
    #create folder
    output_path = os.path.abspath(output_folder)
    output.create_folder(output_path)
    
    if (debug):
        print ("## DEBUG:")
        print ("base:" , base)
        print ("ext:" , ext)
        print ()
        
    gff_parser_caller(gff_file, ref_file, output_path, debug)
    
    
#####
def gff_parser_caller(gff_file, ref_file, output_path, debug):
    
    prot_file = "%s/proteins.fa" % output_path
    csv_file = "%s/df.csv" % output_path            

    with open (ref_file) as in_handle:
        ref_recs = SeqIO.to_dict(SeqIO.parse(in_handle, "fasta"))
    
    if (debug):
        print ("## DEBUG:")
        print (ref_recs)
        print ()

    with open(prot_file, "w") as out_handle:
        SeqIO.write(protein_recs(gff_file, ref_recs, output_path, debug=debug), out_handle, "fasta")
    return (prot_file, csv_file)

#####       
def protein_recs(gff_file, ref_recs, output_path, debug=False):
    #create an empty dataframe. The most important info as first columns
    #all the other file info will be filled next
    columns = ['rec_id', 'locus_tag', 'protein_id', 'gene',  'start',
               'end', 'strand', 'pseudo', 'product', 'Dbxref', 'inference']
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
            if (debug):
                print ("## DEBUG: rec")
                print (rec)
                print ()
        
            for feature in rec.features:
                ## JF
                #if (debug):
                    #print ("## DEBUG: feature")
#                     print (feature)
                    #print ()
# 
#                     ## check feature qualifiers: some might not be always present: pseudo, gene, Name
#                     print (feature.location)
#                     print (feature.location.nofuzzy_start)
#                     print (feature.location.nofuzzy_end)
#                     print (feature.strand)
#                     print (feature.id)
#                     print (feature.qualifiers)
#                     print(feature.type)

                    #print (feature.qualifiers["gene"][0])
                    #print (feature.qualifiers["Name"][0])
                    #print (feature.qualifiers["Parent"][0])
                    #print (feature.qualifiers["locus_tag"][0])
                
                if feature.strand == -1:
                    strand = "neg"
                else:   
                    strand = "pos"
                    
                #we create an ID for each entry     
                protID = feature.type + "_" + rec.id + "_" + str(feature.location.nofuzzy_start) + "_" + str(feature.location.nofuzzy_end) + "_" + strand
#                 if (debug):
#                     print("#DEBUG: protID")
#                     print(protID)
#                     print()
                
                annot_df.loc[protID, ["rec_id", "type", "start", "end", "strand"]] = [ID, feature.type, feature.location.nofuzzy_start, feature.location.nofuzzy_end, strand]
                qualif = feature.qualifiers
#                 if (debug):
#                     print("#DEBUG: Qualifiers")
#                     print(qualif)
#                     print()
                
                for keys, values in qualif.items():
                    
                    #fill the dataframe info
                    if keys == "Note":
                        continue
                    annot_df.loc[protID,[keys]] = [values[0]]
                    
                ## get gene sequence
                gene_seq = Seq.Seq(str(rec.seq[feature.location.nofuzzy_start:feature.location.nofuzzy_end]))

                ## JF
                #if (debug):
                    #print ("## DEBUG: gene_seq")
#                     print (gene_seq)
                    #print()
                        
                
                if feature.type == "CDS":
                    if feature.strand == -1:
                        gene_seq = gene_seq.reverse_complement()
                    
                    #delete STOP symbols
                    protein_seq = gene_seq.translate(to_stop=True)
                    
                    yield(SeqRecord(protein_seq, protID, "", ""))

    csv_file = "%s/df.csv" % output_path            
    annot_df.to_csv(csv_file, header=True)
    csv_length = "%s/length.csv" % (output_path)            
    genome_length.to_csv(csv_length, header=True)
    
    if (debug):
        print("##DEBUG: Dataframe")
        print(annot_df)
        print()
    
    #get genome length for BioCircos plotting  
    genome_length = pd.DataFrame(data=None, columns=["length"])
    ID = rec.id
    length = len(rec.seq)
    genome_length.loc[ID,["length"]]=[length]
    csv_length = "%s/%s_length.csv" % (output_path, rec.id)            
    genome_length.to_csv(csv_length, header=True)
    return(csv_length)
    
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
        
        
        
