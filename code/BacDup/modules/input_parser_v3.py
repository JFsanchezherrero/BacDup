
'''
Created on 16 mar. 2021
@author: alba
'''
############
    # 1: identificar annot_file: GFF / Genbank / None
    
    # 2: Si es genbank: gbf_parser
    
    # 3: Si es GFF:
    # 3.1: Comprobar ref_file: Yes/No
    # 3.2: gff_parser(annot_file, ref_file)
    
    # ToDO: create gtf parser
#############

## useful imports
import os
import sys
import argparse
from Bio import SeqIO

import HCGB
from termcolor import colored

## my modules
import BacDup.scripts.gbf_parser
import BacDup.scripts.gff_parser
import BacDup.scripts.format_checker

##########################
def help_input():
    print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))

##########################
def run_input(arg_dict):
    
    """
    Main function of the input_parser module in BacDup package.
    
    This module prepares data for later gene duplication analysis. 
    
    It allows the user to provide either a single sample, multiple samples or NCBI 
    GenBank IDs to retrieve and obtain the data.    
    """
    
    ## help message
    if arg_dict.help_input():
        help_input()
        exit()
            
    compt = {}
    compt["fasta"] = [".fa", ".faa", ".mpfa", ".fna", ".fsa", ".fas", ".fasta"]
    compt["genbank"] = [".genbank", ".gb", ".gbf", ".gbff", ".gbk"]
    compt["GFF"] = [".gff"]
    
    if arg_dict["debug"]:
        print("## DEBUG: Format input file ")
        print(compt)
        
    ## output
    output_path= os.path.abspath(arg_dict["out_folder"])
    HCGB.functions.file_functions.create_folder(output_path)

    #get  annot file absolute path    
    file_name_abs_path = os.path.abspath(arg_dict["annot_file"])
    #name_file, extension = os.path.splitext(file_name_abs_path)
    
    if arg_dict["debug"]:
        print("## DEBUG: name_file and extension ")
        print(os.path.splitext(file_name_abs_path))
         
    #get protein fasta file and annotation csv file    
    if format_checker.is_gbk(file_name_abs_path) == True:
        if not arg_dict["ref_file"] is None:
            print("######")
            print("Sorry, we don't need a ref file to parse a genbank file")
            print("######")
        else:
        ## call gbf_parser
            prot_file, csv_file = protein_gbf.gbf_parser_caller(arg_dict["annot_file"], output_path, arg_dict["debug"])
            return (prot_file, csv_file)
        if arg_dict["debug"]:
                print("## DEBUG: prot_file ")
                print(prot_file)
    
    elif format_checker.is_gff(file_name_abs_path) == True:
            
        if arg_dict["ref_file"]==None:
            print("######")
            print("Please provide a ref_file FASTA format")
            print("######")
            #print(parser.print_help())
        else:
            ref_name_abs_path = os.path.abspath(arg_dict["ref_file"])
            #ref_name, ref_extension = os.path.splitext(ref_name_abs_path)
            if format_checker.is_fasta(ref_name_abs_path) == True:
                #protein_gff.protein_recs(arg_dict["annot_file"], arg_dict["ref_file"])
                prot_file, csv_file = protein_gff.gff_parser_caller(arg_dict["annot_file"], arg_dict["ref_file"], output_path, arg_dict["debug"])
                return (prot_file, csv_file)
            else:
                print("######")
                print("Compatible ref_file formats: ")
                print(compt["fasta"])
                print("######")  
                 
        if arg_dict["debug"]:
            print("## DEBUG: ref_name and extension ")
            print(os.path.splitext(ref_name_abs_path))
            
    else:
        print("######")
        print("Compatible annot_file formats: ")
        print("#GenBank: ")
        print(compt["genbank"])
        print("#GFF3: ")
        print(compt["GFF"])
        print("######")
       

####################
## arguments
####################
    
## Cuando no se importe como un modulo este script, se ejecutara esto
if __name__ == '__main__':
    parser = ArgumentParser(prog='inputParser',
                            formatter_class=argparse.RawDescriptionHelpFormatter,
                            description="Get proteins sequences from annotation file")
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-a", "--annot_file", metavar="", help="Annotation file: genbank or GFF", required=True)
    requiredNamed.add_argument("-o", "--out_folder", metavar= "", help="Results folder", required=True)
    parser.add_argument("-r", "--ref_file", metavar="", help="Genome references FASTA file")
    #parser.add_argument("-p", "--prot_file", metavar="", help="Protein sequence file")
    parser.add_argument("--debug", action="store_true", default=False)   
    
    arg = parser.parse_args()
    arg_dict = vars(arg)
    
  
    if arg.annot_file is None or arg.out_folder is None:
        print("######")
        print("Please provide either an annotation file or/and an output path folder")
        print("######")
        print(parser.print_help())
        exit()
        
    if arg.debug:
        print("##DEBUG: ##")
        print("arguments dictionary: ")
        print(arg)

    input_parser(arg_dict)

if __name__ == '__main__':
    pass