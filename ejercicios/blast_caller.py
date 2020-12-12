'''
Created on 28 nov. 2020

@author: alba
'''

import os
import sys

import subprocess
import argparse
from argparse import ArgumentParser

#from Bio.Blast.Applications import NcbiblastpCommandline

#https://github.com/HCGB-IGTP/HCGB_python_functions/blob/ce1762b709e3b9094c0630f63bfb9d33544ed4c3/HCGB/functions/blast_functions.py
#####
def makeblastdb(makeblastdb_exe, fasta_file, db_path_name):
    print()
    
    cmd_makeblastdb = "%s -in %s -input_type fasta -dbtype %s -out %s" %(makeblastdb_exe, fasta_file, 'prot', db_path_name)
    code = system_call(cmd_makeblastdb)
        ##
        ## -> ¿distinguir que sea un fasta con aminoácidos y no nucleótidos?
        ##
    if (code == 'FAIL'):
        print ('****ERROR: Some error happened during the makeblastDB command')
        print (cmd_makeblastdb)
        exit()
    else:
        return(True)

#####        
def blastp_caller(blastp_exe, fasta_file, db_path_name, output_file):
    print()
    cmd_makeblast_results = "%s -query %s -db %s -outfmt \'6 std qlen slen\' -num_threads 1 -out %s" %(blastp_exe, fasta_file, db_path_name, output_file)
    code_results = system_call(cmd_makeblast_results)
    
    if (code_results == 'FAIL'):
        print ('****ERROR: Some error happened during the blastp command')
        print (cmd_makeblast_results)
        exit()
    else:
        return(True)
    
#####    
def create_blast_results(arg_dict):
    
    #possible extensions    
    compt = {}
    compt["fasta"] = [".fa", ".faa", ".mpfa", ".fna", ".fsa", ".fas", ".fasta"]
    
    ## generate blastdb for genome
    #phr is the header file, pin is the index file, psq is the sequence file
    file_name_abs_path = os.path.abspath(arg_dict["fasta_file"])
    name_file, extension = os.path.splitext(file_name_abs_path)
    basename= os.path.basename(name_file)
    #output_path = "%s_blastp_results.txt" % arg_dict["db_name"]
    if arg.debug:
        print("## Debug: name_file and extension ")
        print(os.path.splitext(file_name_abs_path))
            
    if extension in compt["fasta"]:
            
        if arg.debug:
            print("## DEBUG: Format input file ")
            print(compt)
            
        ## create blast database
        makeblastdb_exe = arg_dict["blast_folder"] + "/makeblastdb"
        blastp_exe = arg_dict["blast_folder"] + "/blastp"
        ##-> blast_folder = /usr/bin
        #name_file = "name"
        if (arg_dict["db_name"]):
            db_path_name = os.path.abspath(arg_dict["db_name"])+ "/" + arg_dict["db_name"]
            output_file = os.path.abspath(arg_dict["db_name"]) + "/BLAST_raw_results.txt"
        elif (arg_dict["out_folder"]):
            db_path_name = os.path.abspath(arg_dict["out_folder"]) + "/" + basename + "_db"
            output_file = os.path.abspath(arg_dict["out_folder"]) + "/" + basename + "_BLAST_raw_results.txt"
        else:
            db_path_name = basename + "_db"
            output_file = "BLAST_raw_results.txt"
    
        if (os.path.isfile(db_path_name + '.phr')):
            print ("+ BLAST database is already generated...")
        else:
            makeblastdb(makeblastdb_exe, arg_dict["fasta_file"], db_path_name)
        
        ## create blastp outfile
            
            blastp_caller(blastp_exe, arg_dict["fasta_file"], db_path_name, output_file)
                 
    else:
        print("#####")
        print("Please provide a FASTA file")
        print(compt)
        print("#####")
            
        

#https://github.com/HCGB-IGTP/HCGB_python_functions/blob/2b3b1132fb885c8cb22f4d10fd4c00c25fa10fb8/HCGB/functions/system_call_functions.py
#####
def system_call(cmd, returned=False, message=True):
    """Generates system call using subprocess.check_output"""
    ## call system
    ## send command
    if (message):
        print ("[** System: %s **]" % cmd)

    try:
        out = subprocess.check_output(cmd, shell = True)
        if (returned):
            return (out)
        return ('OK')
    except subprocess.CalledProcessError as err:
        if (returned):
            return (err.output)
        if (message):
            print ("** ERROR **")
            print (err.output)
            print ("** ERROR **")
        
        return ('FAIL')
    

    
####################
## Arguments
####################
    

parser = ArgumentParser(prog='makeblastDB',
                        formatter_class=argparse.RawDescriptionHelpFormatter,
                        description="Create a BLAST database")
parser.add_argument("-d", "--db_name", metavar="", help="New database name")
parser.add_argument("-f", "--fasta_file", metavar="", help="Protein sequences FASTA file")
parser.add_argument("-b", "--blast_folder", metavar="", help="BLAST binary folder")
# parser.add_argument("-ev", "--evalue", type=float, metavar="", help="BLAST e-value: number of expected hits of similar quality (score) that could be found just by chance.")
# parser.add_argument("-b", "--bitscore", type=int, metavar="", help="BLAST bit-score: requires size of a sequence database in which the current match could be found just by chance.")
parser.add_argument("-o", "--out_folder", metavar= "", help="Results folder")
parser.add_argument("--debug", action="store_true", default=False)   

arg = parser.parse_args()
arg_dict = vars(arg)
    
if arg.fasta_file is None:
    print("#####")
    print("Please provide a proteins sequences FASTA file")
    print("#####")
    print(parser.print_help())
     
if arg.blast_folder is None:
    print("#####")
    print("Please provide BLAST binary folder")
    print("#####")
    print(parser.print_help())
         
if arg.debug:
    print("##DEBUG: ##")
    print("arguments dictionary: ")
    print(arg)
  
create_blast_results(arg_dict) 
  
    