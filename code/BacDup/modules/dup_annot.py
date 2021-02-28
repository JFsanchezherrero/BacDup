'''
Created on 13 dic. 2020

@author: alba
'''

import os
import argparse
from argparse import ArgumentParser


import pandas as pd
from collections import defaultdict 


import dup_searcher
import input_parser


#############
### DONE! ###
#############
# 1. provide a GFF+FASTA/GBF file
#
# 2. get sequence protein FASTA file ############
#                                    ############## input_parser --> protein_gbf or protein_gff
#    get annotation csv file from GFF/GBF file ##
#
# 3. convert protein FASTA file into a BLAST database ##
# 4. get duplicated proteins BLASTA .txt file ############ dup_searcher --> blast_caller
# 5. get filtered results into a .csv file #############
###############
#### TO DO ####
###############
# 6. get results annotation info from annotation csv file ######
#                                          ######################## dup_annot -> dup_searcher -> input_parser  
# 7. classify results by duplicated number #####################
###################################################################################

# def get_seq:
#     # get seq proteins FASTA file
#     # get annot csv file
#     input_parser.input_parser(arg_dict)
#     
# def get_dup:
    
# def get_seq:
#     # get seq proteins FASTA file
#     # get annot csv file
#     input_parser.input_parser(arg_dict)
#     
# def get_dup:
def get_dupannot(arg_dict):
    
    '''
    Get an annotation file for duplicated proteins
    '''
    
    compt = {}
    compt["fasta"] = [".fa", ".faa", ".mpfa", ".fna", ".fsa", ".fas", ".fasta"]
    compt["genbank"] = [".genbank", ".gb", ".gbf", ".gbff", ".gbk"]
    compt["GFF"] = [".gff"]
#     
    output_path = os.path.abspath(arg_dict["out_folder"]) # ERROR: check user provides output
#     
#   #user provides an GFF/GBK file
    if arg_dict["text_file"] is None:
        if arg_dict["annot_file"]:
            file_name_abs_path = os.path.abspath(arg_dict["annot_file"])
            name_file, extension = os.path.splitext(file_name_abs_path)
            if extension in compt["GFF"]:
                if arg_dict["ref_file"] is None:
                    print("######")
                    print("Please provide a ref_file FASTA format")
                    print("######")
                    print(parser.print_help())
                    exit()
            prot_file, csv_file = input_parser.input_parser(arg_dict)
            arg_dict["fasta_file"] = prot_file
            arg_dict["annot_table"] = csv_file
            
        elif arg_dict["fasta_file"]:
            if arg_dict["annot_table"] is None:
                print("#####")
                print("Please provide an annotation table")
                print("#####")
                print(parser.print_help())
                exit()    
            else:
                file_name_abs_path = os.path.abspath(arg_dict["annot_table"])
                name_file, extension = os.path.splitext(file_name_abs_path)
                if extension != ".csv":
                    print("#####")
                    print("Please provide a .csv file")
                    print("#####")
                    print(parser.print_help())
                    exit()
            ## ToDo check annot_table and fasta_file headers are the same ##

        # ERROR: check user provides BLAST folder if required

## create blast results
        filtered_data = dup_searcher.filter_data(arg_dict)
        annot_table = pd.read_csv(csv_file, index_col=0)
    
    else:
        # get filtered results file
        filtered_data = dup_searcher.filter_data(arg_dict)

        # ERROR: check user provides Annotation CSV

        annot_table = pd.read_csv(arg_dict["annot_table"], index_col=0)
        ## debug 
        #print (annot_table)
    #get duplicated protein list
    qseqid = list(filtered_data["qseqid"])
    sseqid =list(filtered_data["sseqid"])
    qseqid.extend(sseqid)
    prot_id = list(set(qseqid))
    
    #get filtered_annot table
    filtered_annot = annot_table.loc[prot_id]

    # keep or maintain pseudogenes
    if (arg_dict['pseudo']):
        # use them
        print ("+ Pseudogenes would be used in the analysis")
    else:
        filtered_annot = filtered_annot.loc[filtered_annot['pseudo'] != True] 

    ## debug 
    #print (filtered_annot)    
    dup_annot = get_dup(filtered_data, filtered_annot)

    dup_annot_file = "%s/dup_annot.csv" % output_path
    dup_annot.to_csv(dup_annot_file, header=True)
    print("#####")
    print("Found %s groups of duplicates with a total of %s proteins duplicated from %s proteins on the original file" % (len(dup_annot["dup_id"].unique()), dup_annot.shape[0], annot_table.shape[0]))
    print("#####")
    return(dup_annot)


def get_dup(blast_results_df, dup_annot_df):
    
    ''' 
    get duplicated id for each duplicated protein on annotation table
    
    '''
    
    # 1st round
    relations_dict = defaultdict(list) 
    for index, row in blast_results_df.iterrows():
        relations_dict[row['qseqid']].append(row['sseqid'])
    
    ## debug
    # print (relations_dict)

    ## 2nd round
    new_relations_dict = defaultdict(list)
    dups=0
    for key, value in relations_dict.items():
        ## debug
        ## print ()
        ## debug
        ## print ("key: " + key)
        ## debug
        ## print ("value: " + str(value))

        stop=False
        for dup_id, new_value in new_relations_dict.items():
            if key in new_value:
                stop=True
                ## debug
                ## print ("Belongs to group: " + dup_id)

        if not stop:
            for key2, value2 in relations_dict.items():
                if (key == key2):
                    continue
                else:
                    if (key2 in value): 
                        for i in value2: 
                            if i not in value: 
                                value.extend(i)

            dups += 1
            value.append(key)
            new_relations_dict[str(dups)] = value
            ## debug
            ## print(new_relations_dict)
            ## debug
            ## print("**")
            
    ## print ()
    ## print (new_relations_dict)
    ## print ()

    ## Create data
    df_data = pd.DataFrame(columns=('index', 'dup_id'))
    for dup_id, new_value in new_relations_dict.items():
        for i in new_value:
            df_data.loc[i] = (i, dup_id)

    ## merge information    
    dup_annot_df = dup_annot_df.join(df_data)
    dup_annot_df = dup_annot_df.drop(columns='index')
    return (dup_annot_df)            
        

############    
parser = ArgumentParser(prog='dupAnnotation',
                        formatter_class=argparse.RawDescriptionHelpFormatter,
                        description="Get an annotation file with duplicated protein on genome.")

annotparser = parser.add_argument_group('Annotation parser options named arguments')
annotparser.add_argument("-a", "--annot_file", metavar="", help="Annotation file: genbank or GFF.")
annotparser.add_argument("-r", "--ref_file", metavar="", help="Genome references FASTA file.")
###
blastoptions = parser.add_argument_group('BLAST options named arguments')
blastoptions = parser.add_argument_group('BLAST options named arguments')
blastoptions.add_argument("-b", "--blast_folder", metavar="", default="/usr/bin", help="BLAST binary folder. **Note blast_folder=/usr/bin is set by default**")    
blastoptions.add_argument("-bs", "--bitscore", type=float, metavar="", default=50, help="BLAST bit-score: requires size of a sequence database in which the current match could be found just by chance. **Note bit_score = 50 is set by default**")
blastoptions.add_argument("-d", "--db_name", metavar="", help="New database name")
blastoptions.add_argument("-e", "--evalue", type=float, metavar="", default= 1e-05, help="BLAST e-value: number of expected hits of similar quality (score) that could be found just by chance. **Note e-value = 1e-05 is set by default**")
blastoptions.add_argument("-p", "--percentage", type=float, metavar="", default=85, help="Percentage of alignment in query. *Note pident = 85 is set by default**")
blastoptions.add_argument("-pi", "--pident", type=int, metavar="", default=85, help="Percentage of similarity in alignment. **Note percentage = 85 is set by default**")

parser.add_argument("-c", "--annot_table", metavar="", help="Genome annotation .csv file previously analyzed.")

parser.add_argument("-t", "--text_file", metavar="", help="Blast raw results text file.")
###
parser.add_argument("--pseudo", action="store_true", default=False, help="Wether to use pseudogenes or not")
###
parser.add_argument("-o", "--out_folder", metavar= "", help="Results folder")
parser.add_argument("--debug", action="store_true", default=False)

arg = parser.parse_args()
arg_dict = vars(arg)
print("######")
print("e-value = %s" % arg.evalue)
print("bitscore = %s" % arg.bitscore)
print("percentage alignment in query = %s" % arg.percentage)
print("percentage of identity = %s" % arg.pident)
print("blast_folder = %s" % arg.blast_folder)
print("######")

if arg.annot_file is None and arg.fasta_file is None and arg.text_file is None:
    print("######")
    print(parser.print_help())
    print("######")
    exit() 

if arg.debug:
    print("##DEBUG: ##")
    print("arguments dictionary: ")
    print(arg)


if __name__ == '__main__':
    get_dupannot(arg_dict)
    


