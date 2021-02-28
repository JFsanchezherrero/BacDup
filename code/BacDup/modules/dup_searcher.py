'''
Created on 3 dic. 2020

@author: alba
'''
import os

import pandas as pd
import numpy as np

import argparse
from argparse import ArgumentParser

import blast_caller

#https://biopython-tutorial.readthedocs.io/en/latest/notebooks/07%20-%20Blast.html#Parsing-BLAST-output

## -> evalue col 11
## bitscore -> col 12

  
#test_file = "/home/alba/git/TFM_UOC_AMoya/test_BLAST_raw_results.txt"  
######
#test_file:
#Total:16prots
#A=B=C; D=E=F (34aa changed);
#G=H=I (180aa from A + whole J);
#J=K=L; M=N=O (17aa changed);
#P

######
def filter_data(arg_dict):
    
    '''
    get duplicated proteins and generate a filtered and sort dataframe depending on arguments values
    
    if a BLAST results text file is provided, it should have next columns names:
    
    #################################################################
    # qseqid -> query (e.g., unknown gene) sequence id
    # sseqid -> subject (e.g., reference genome) sequence id
    # pident -> percentage of identical matches
    # length -> alignment length (sequence overlap)
    # mismatch -> number of mismatches
    # gapopen -> number of gap openings
    # qstart -> start of alignment in query
    # qend -> end of alignment in query
    # sstart -> start of alignment in subject
    # send -> end of alignment in subject
    # evalue -> expect value
    #            number of expected hits of similar quality (score) that could be found just by chance.
    #            Blast results are sorted by E-value by default (best hit in first line).
    # bitscore -> bit score
    #            requires size of a sequence database in which the current match could be found just by chance.
    #            The higher the bit-score, the better the sequence similarity
    # qlen -> Query sequence length
    # slen -> Subject sequence length
    # aln_pct_qlen -> alignment percentage in query (length/qlen*100)
    # aln_pct_slen -> alignment percentage in subject (length/slen*100)
    ################################################################
    
    '''
    columns = ["qseqid", "sseqid", "pident", "length",
           "mismatch", "gapopen", "qstart", "qend",
           "sstart", "send", "evalue", "bitscore",
           "qlen", "slen", "aln_pct_qlen", "aln_pct_slen"]

    if arg_dict["text_file"]:
        if arg_dict["fasta_file"]==None:
            
            file_name_abs_path = os.path.abspath(arg_dict["text_file"])
            name_file, extension = os.path.splitext(file_name_abs_path)
            
            if extension == ".txt":
                #dataframe with raw blast results
                raw_blast = pd.read_csv(file_name_abs_path, sep="\t",
                                        header = None, names=columns)
                
		## ERROR: remove to allow users to filter again

		#if (arg_dict["out_folder"]):
                #    output_path= os.path.abspath(arg_dict["out_folder"])
                #    sort_csv = "%s/filtered_results.csv" % output_path
                #else:
                #    sort_csv = "filtered_results.csv"
                #return (sort_csv)    

            ## TODO verify columns
            
            else:
                print("#####")
                print("Please provide a BLAST results .txt file")
                print("#####")
                exit()
        else:
            print("#####")
            print("Please provide either a BLAST results .txt file OR an annotation FASTA file")
            print("#####")
            exit()
            
    elif arg_dict["fasta_file"]:
        if arg_dict["blast_folder"] is None:
            print("#####")
            print("Please provide BLAST binary folder")
            print("#####")
            #print(parser.print_help())
            exit()
        
        else:    
            raw_blast = blast_caller.create_blast_results(arg_dict)
            raw_blast = pd.read_csv(raw_blast, sep="\t", header = None, names=columns)
            
#fill aln_pct columns
    raw_blast["aln_pct_qlen"] = (raw_blast["length"]/raw_blast["qlen"]*100).round(2)
    raw_blast["aln_pct_slen"] = (raw_blast["length"]/raw_blast["slen"]*100).round(2)

#filter mirror pairs
    del_mirror =pd.DataFrame(np.sort(raw_blast[["qseqid", "sseqid"]], axis=1))
    raw_blast = raw_blast[~del_mirror.duplicated()]
          
#when Blast_file_results is done, evalue>10 is gone
    filter_pident= raw_blast["pident"] >= arg_dict["pident"]

    filter_evalue = raw_blast["evalue"] <= arg_dict["evalue"]
    filter_bitscore = raw_blast["bitscore"] >= arg_dict["bitscore"]
    filter_alnpct_slen = raw_blast["aln_pct_slen"] >= arg_dict["percentage"]
    filter_alnpct_qlen = raw_blast["aln_pct_qlen"] >= arg_dict["percentage"]
    filter_id = raw_blast["qseqid"] != raw_blast["sseqid"]

    filtered_results = raw_blast[filter_pident & filter_evalue & filter_bitscore & filter_alnpct_qlen & filter_alnpct_slen & filter_id]
    
#sort by aln_pct_qlen (desc), evalue(asc), bitscore(desc)
    by_alnpct = filtered_results.sort_values(["aln_pct_qlen", "evalue", "bitscore"],
                                   ascending=[False, True, False])


#     print(by_alnpct)
#     print("#####")
#     print("%s pairs filtered from %s" % (by_alnpct.shape[0], raw_blast.shape[0]))
#     print("#####")

#save results as a .csv file
    if (arg_dict["db_name"]): 
        output_path= os.path.abspath(arg_dict["db_name"])
        sort_csv = "%s/filtered_results.csv" % output_path
    elif (arg_dict["out_folder"]):
        output_path= os.path.abspath(arg_dict["out_folder"])
        sort_csv = "%s/filtered_results.csv" % output_path
    else:
        sort_csv = "filtered_results.csv"
    
    by_alnpct.to_csv(sort_csv, header=True, index=False)
    return(by_alnpct)     
    
    
######
##
######

if __name__ == '__main__':
    parser = ArgumentParser(prog='dupSearcher',
                            formatter_class=argparse.RawDescriptionHelpFormatter,
                            description="Search for genomic duplicated proteins")
    parser.add_argument("-t", "--text_file", metavar="", help="Blast raw results text file")
    parser.add_argument("-f", "--fasta_file", metavar="", help="Protein sequences FASTA file")
    parser.add_argument("-o", "--out_folder", metavar= "", help="Results folder")
    parser.add_argument("--debug", action="store_true", default=False)
    ##
    blastoptions = parser.add_argument_group('BLAST options named arguments')
    blastoptions.add_argument("-b", "--blast_folder", metavar="", help="BLAST binary folder. **Note blast_folder=/usr/bin is set by default**")    
    blastoptions.add_argument("-bs", "--bitscore", type=float, metavar="", default=50, help="BLAST bit-score: requires size of a sequence database in which the current match could be found just by chance. **Note bit_score = 50 is set by default**")
    blastoptions.add_argument("-d", "--db_name", metavar="", help="New database name")
    blastoptions.add_argument("-e", "--evalue", type=float, metavar="", default= 1e-05, help="BLAST e-value: number of expected hits of similar quality (score) that could be found just by chance. **Note e-value = 1e-05 is set by default**")
    blastoptions.add_argument("-p", "--percentage", type=float, metavar="", default=85, help="Percentage of alignment in query. *Note pident = 85 is set by default**")
    blastoptions.add_argument("-pi", "--pident", type=int, metavar="", default=85, help="Percentage of similarity in alignment. **Note percentage = 85 is set by default**")
    
    arg = parser.parse_args()
    arg_dict = vars(arg)
    filter_data(arg_dict)
 

     
    if arg.evalue is None:
        print("#####")
        print("Note e-value = 1e-05 is set by default")
        print("#####")
        print(parser.print_help())
         
    if arg.bitscore is None:
        print("#####")
        print("Note bit_score = 50 is set by default")
        print("#####")
        print(parser.print_help())
        
    if arg.percentage is None:
        print("#####")
        print("Note alignment in query percentage = 80% is set by default")
        print("#####")
        print(parser.print_help())
    
    if arg.debug:
        print("##DEBUG: ##")
        print("arguments dictionary: ")
        print(arg)
