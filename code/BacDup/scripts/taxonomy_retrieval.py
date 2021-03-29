#!/usr/bin/env python3
"""
Created by Joe R. J. Healey; Nick Youngblut
Modified in March 2021 by Jose F. Sanchez-Herrero

Original code retrieved from ncbi-genome-download/contrib/gimme_taxa.py

Perform various queries and tasks with the NCBI taxonomy database,
via the ETE3 toolkit.

"""

## https://academic.oup.com/database/article/doi/10.1093/database/baaa062/5881509
## https://www.ncbi.nlm.nih.gov/books/NBK179288/
## https://github.com/kblin/ncbi-genome-download
## https://github.com/kblin/ncbi-genome-download/blob/master/contrib/gimme_taxa.py

import sys
import os
import argparse
from ete3 import NCBITaxa
from HCGB.functions.aesthetics_functions import debug_message
from HCGB.functions import time_functions
import HCGB
from termcolor import colored
import random

## conda install -c etetoolkit ete3 ete_toolchain

######################################################################
def parse_taxid(tax_id, ncbi, option, debug):
    if (option=="info"):
        print()
        (tax_name, taxid, rank, lineage) = taxon_info(tax_id, ncbi, debug)
        
        print ("Result:")
        print ("Name: " + tax_name)
        print ("Rank: " + rank)
        print ("Taxid: " + str(taxid))
        list_lineage = lineage.split(";")
        print ("Lineage:")
        
        for tax in list_lineage:
            tax_split = tax.split(":")
            print ("\t" + '{}\t{}'.format(tax_split[0], tax_split[1])) 
        
        ## return info
        return (tax_name, taxid, rank, lineage)
        
        
    elif (option=="unravel"):
        info = unravel_taxid(tax_id, ncbi, debug)
        return info

######################################################################
def unravel_taxid(tax_id, ncbi, debug):
    
    ## check the rank provided    
    (tax_name, taxid, rank, lineage) = taxon_info(tax_id, ncbi, debug)
    
    if debug:
        debug_message("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        debug_message('tax_name: '+ tax_name, "yellow")
        debug_message('taxid:' + str(taxid) , "yellow")
        debug_message('rank: ' + rank, "yellow")
        debug_message('lineage: ' + lineage, "yellow")
    
    ##
    list_taxids = []
    list_lineage = lineage.split(";")
    for tax in list_lineage:
        tax_split = tax.split(":")
        ## check the rank provided    
        (tax_name2, taxid2, rank2, lineage2) = taxon_info(tax_split[0], ncbi, debug)
        if (rank2 == "species"):
            list_taxids.append(taxid2)
        elif (rank2 == "serotype"):
            list_taxids.append(taxid2)
            
    return (list_taxids)   

        
######################################################################
def init_db_object(debug):
    """Instantiate the ete taxonomy object     
    Created by Joe R. J. Healey; Nick Youngblut
    Original code.
    """
    # Instantiate the ete NCBI taxa object
    print ("+ ------------------------------------- +")
    print ("+ Looking for NCBI taxonomy database:")
    ncbi = NCBITaxa()
    
    ## dbfile location
    if debug:
        debug_message("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        debug_message('NCBI Taxonomy database is stored at {}\n'.format(ncbi.dbfile), "yellow")
    
    ## folder would be download here: ~/.etetoolkit/taxa.sqlite 
    db_folder = os.path.dirname(format(ncbi.dbfile))
    
    ## check timestamp, update if necessary
    filename_stamp_parse = db_folder + '/timestamp_db.txt'
    if os.path.isfile(filename_stamp_parse):
        stamp = time_functions.read_time_stamp(filename_stamp_parse)
        days_passed = time_functions.get_diff_time(stamp)
        
        ## debug messages
        if debug:
            debug_message('Database previously initiated', "yellow")
            debug_message('on date: {}'.format(stamp), "yellow")
            debug_message('Days passed: {}'.format(days_passed), "yellow")
        
        if (days_passed > 30):
            ## update_db
            update_db(ncbi, db_folder, debug)
        else:
            ## debug messages
            if debug:
                debug_message('No need to update db', "yellow")
                debug_message("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

            print (colored("\tA previous command generated results on: %s [%s]" %(stamp, 'init database'), 'yellow'))
    else:
        ## create first timestamp
        time_functions.print_time_stamp(filename_stamp_parse)

    return ncbi

######################################################################
def update_db(ncbi_db, db_folder, debug):
    """Update database
    Created by Joe R. J. Healey; Nick Youngblut
    Original code.
    """
    
    ## debug messages
    if debug:
        debug_message('Update database at {}\n'.format(ncbi_db.dbfile), "yellow")
            
    print ('Updating the taxonomy database. This may take several minutes...\n')
    ncbi_db.update_taxonomy_database()
    
    ## print timestamp
    filename_stamp_parse = os.path.abspath(db_folder + '/timestamp_db.txt')
    time_functions.print_time_stamp(filename_stamp_parse)

    return ncbi_db

######################################################################
def desc_taxa(taxid, ncbi, outFH):
    """Write descendent taxa for taxid
    Created by Joe R. J. Healey; Nick Youngblut
    Slightly modified. 
    
    Returns python dictionary with descendant taxid and name.
    """
    # Main feature of the script is to get all taxa within a given group.    
    descendent_taxa = ncbi.get_descendant_taxa(taxid)
    descendent_taxa_names = ncbi.translate_to_names(descendent_taxa)

    dict_Descent = {}
    for dtn, dt in zip(descendent_taxa_names, descendent_taxa):
        dict_Descent[dt] = dtn
    
    return dict_Descent

######################################################################
def taxon_info(taxid, ncbi, debug):
    """Write info on taxid
    Originally created by Joe R. J. Healey; Nick Youngblut
    
    Slightly modified code by JFSanhcezherrero
    """
    try:
        taxid = int(taxid)
        tax_name = ncbi.get_taxid_translator([taxid])[taxid]
        rank = list(ncbi.get_rank([taxid]).values())[0]
        lineage = ncbi.get_taxid_translator(ncbi.get_lineage(taxid))
        lineage = ['{}:{}'.format(k,v) for k,v in lineage.items()]
        lineage = ';'.join(lineage)
        return(tax_name, taxid, rank, lineage)
    except:
        print (colored("** ERROR: ", 'red'))
        print ("No valid tax id provided: " + str(taxid))
        exit()
        return("", "", "", "")

######################################################################
def name2taxid(taxids, ncbi):
    """Converting taxon names to taxids
    Created by Joe R. J. Healey; Nick Youngblut
    Original code.
    """
    
    # If names were provided in taxid list, convert to taxids
    #args.taxid = args.taxid.replace('"', '').replace("'", '').split(',')
    #args.taxid = name2taxid(args.taxid, ncbi)
    
    new_taxids = []
    for taxid in taxids:
        try:
            new_taxids.append(ncbi.get_name_translator([taxid])[taxid][0])
        except KeyError:
            try:
                new_taxids.append(int(taxid))
            except ValueError:
                msg = 'Error: cannot convert to taxid: {}'
                raise ValueError(msg.format(taxid))

    return new_taxids

######################################################################
def help_options():
    print ("\nUSAGE: python %s db_folder taxID info/unravel\n"  %os.path.realpath(__file__))


######################################################################
def main():

    ## control if options provided or help
    if len(sys.argv) > 1:
        print ("")
    else:
        help_options()
        exit()
    
    ## get arguments provided
    taxid = sys.argv[1]
    option = sys.argv[2]
    debug=True
    
    ## init database
    ncbi = init_db_object(debug)
    
    ## parse
    info = parse_taxid(taxid, ncbi, option, debug)
    
    if debug:
        debug_message("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        debug_message('info\n', "yellow")
        print(info)
        debug_message("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        
    print ("+ ------------------------------------- +\n\n")

    ## header
    # taxon_info: ['name', 'taxid', 'rank', 'lineage']
    # else:
    # ['parent_taxid', 'descendent_taxid','descendent_name']
    
    ## body
    #for taxid in args.taxid:
    #    if args.taxon_info:
    #        taxon_info(taxid, ncbi, outFH)
    #    else:
    #        desc_taxa(taxid, ncbi,  outFH, args.just_taxids)
            
    
############################################################
if __name__== "__main__":
    main()

