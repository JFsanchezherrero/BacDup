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

from BacDup.scripts.NCBI_downloader import NCBI_get_info

## conda install -c etetoolkit ete3 ete_toolchain

######################################################################
def parse_taxid(tax_id, ncbi, option, debug):
    """Function to parse according to option: provide info or return unravelled data    
    """
    ############
    ## debug messages
    ############
    if debug:
        debug_message('parse_taxid:', "yellow")
        if tax_id.isdigit():
            debug_message('tax_id: ' + str(tax_id), "yellow")
        else:
            debug_message('tax_id: ' + tax_id, "yellow")
            debug_message('conversion needed: ', "yellow")
             
        debug_message('option: ' + option, "yellow")
    
    ############
    ## convert to tax ID
    ############
    if not tax_id.isdigit():
        ## convert name to taxid integer
        print ("+ Convert to NCBI taxonomy ID")
        print ("\tSource: " + tax_id)
        tax_id = name2taxid([tax_id], ncbi)
        (tax_name, taxid, rank, lineage) = taxon_info(tax_id, ncbi, debug)
        print ("\tRank: " + rank)
        print ("\tID: " + str(taxid))
        
    ############
    ## parse accordingly
    ############
    if (option=="info"):
        print()
        (tax_name, taxid, rank, lineage) = taxon_info(tax_id, ncbi, debug)
        
        print ("----------------------------------------------")
        print ("Result:")
        print ("Name: " + tax_name)
        print ("Rank: " + rank)
        print ("Taxid: " + str(taxid))
        list_lineage = lineage.split(";")
        print ("Lineage:")
        
        for tax in list_lineage:
            tax_split = tax.split(":")
            print ("\t" + '{}\t{}'.format(tax_split[0], tax_split[1])) 
        
        print ("----------------------------------------------")
        print ()
        
        ## return info
        return (tax_name, taxid, rank, lineage)
        
    ############
    ### call unravel taxid information
    ############
    elif (option=="unravel"):
        return(unravel_taxid(tax_id, ncbi, debug))
        
###########################################################################
def get_GenBank_ids(data_folder, taxID_list, random_k, debug, assembly_level_given='complete', group_given='bacteria', section_given='genbank'):
    '''
    This function retrieves information from GenBank for the taxids provided
    using ncbi_genome_download module in NCBI_downloader script
    
    :param data_folder: Database folder containing entries
    :param taxID_list: List of taxids to retrieve
    :param random_k: Number of entries to randomly select. Provide -1 for all entries.
    :param debug: True/False for debugging messages
    :param assembly_level_given: Either complete, scaffold, chromosome or all.
    :param group_given: Either bacteria, archaea, virus or eukarya
    
    :type data_folder: str
    :type taxID_list: list
    :type random_k: int
    :type debug: bool
    :type assembly_level_given: str
    :type group_given: str
    
    :return dict
    '''

    len_list = len(taxID_list)

    ## debug messages
    if debug:
        debug_message('get_GenBank_ids:', "yellow")
        debug_message('taxID_list: %s items' %len_list, "yellow")
        print (taxID_list)
        debug_message('group: ' + group_given, "yellow")
        debug_message('assembly level: ' + assembly_level_given, "yellow")
        debug_message('selection: ' + str(random_k), "yellow")
        
        
    ## retrieve NCBI GenBanks ids with com
    dict_entries = NCBI_get_info(section_given, data_folder, taxID_list, debug, assembly_level_given, group_given)
    dict_entries_len = len(list(dict_entries.keys()))
    
    ##
    if random_k<0:
        list_entries = list(dict_entries.keys())
        print ('All %s entries selected' %dict_entries_len)
    
    else:
        print ("Selecting random entries retrieved:")
        list_entries = random.choices(list(dict_entries.keys()), k=random_k)
        print ('%s entries selected out of %s' %(str(random_k), str(dict_entries_len)))
    
    ## debug messages
    if debug:        
        debug_message('list_entries:', "yellow")
        print (list_entries)
        debug_message("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        
    ##
    return(list_entries, list(dict_entries.keys()))


######################################################################
def unravel_taxid(tax_id, ncbi, debug):
    """This function unravels information and obtains children taxids for each taxid.
    
    If taxid corresponds to serotype or species, no further processing is done. On the 
    other hand, if genes, family, order or any other rank is provided, all subranks would 
    retrieved. It also takes into account serotypes and accomodates information.
    
    It returns a list of all taxids included within the tax_id provided.
    """
    
    ## check the rank provided    
    (tax_name, taxid, rank, lineage) = taxon_info(tax_id, ncbi, debug)
    
    ## debug messages
    if debug:
        debug_message('tax_name: '+ tax_name, "yellow")
        debug_message('taxid:' + str(taxid) , "yellow")
        debug_message('rank: ' + rank, "yellow")
        debug_message('lineage: ' + lineage, "yellow")
    
    ##
    list_taxids = []
    
    ## taxid provided is either a serotype or strain: directly to retrieve
    if (rank in ("species", "serotype", "strain")):
        list_taxids.append(taxid)
    else:
        
        ## get descendant
        dict_descent = desc_taxa(taxid, ncbi, debug)
        for tax, name in dict_descent.items():
            
            ## add taxa retrieved
            list_taxids.append(tax)
            
            ## check the rank provided and decompose    
            (tax_name2, taxid2, rank2, lineage2) = taxon_info(tax, ncbi, debug)
            list_lineage = lineage2.split(";")
            for tax3 in list_lineage:
                tax_split = tax3.split(":")
                ## check the rank provided: add also species or serotype
                (tax_name3, taxid3, rank3, lineage3) = taxon_info(tax_split[0], ncbi, debug)
                if (rank3 in ("species", "serotype")):
                    list_taxids.append(taxid3)
            
    ## return uniq list of ids
    return (list(set(list_taxids)))
        
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
def desc_taxa(taxid, ncbi, debug):
    """Write descendent taxa for taxid
    Created by Joe R. J. Healey; Nick Youngblut
    Slightly modified. 
    
    Returns python dictionary with descendant taxid and name.
    """
    ## debug messages
    if debug:
        debug_message("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        debug_message('desc_taxa: '+ str(taxid), "yellow")
    
    # Main feature of the script is to get all taxa within a given group.    
    descendent_taxa = ncbi.get_descendant_taxa(taxid)
    descendent_taxa_names = ncbi.translate_to_names(descendent_taxa)

    dict_Descent = {}
    for dtn, dt in zip(descendent_taxa_names, descendent_taxa):
        dict_Descent[dt] = dtn
    
    ## debug messages
    if debug:
        debug_message('dict_Descent: ', "yellow")
        print (dict_Descent)
        debug_message("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    
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


######################################################################
def get_superKingdom(tax_id, ncbi, debug):
    """For a given tax_id get superkingdom from NCBI taxonomy ID
    """
    
    if debug:
        debug_message("get_superKingdom: ", 'yellow')
        debug_message("tax_id: " + str(tax_id), 'yellow')
    
    (tax_name2, taxid2, rank2, lineage2) = taxon_info(tax_id, ncbi, debug)
    list_lineage = lineage2.split(";")
    for tax3 in list_lineage:
        tax_split = tax3.split(":")
        ## check the rank provided: add also species or serotype
        (tax_name3, taxid3, rank3, lineage3) = taxon_info(tax_split[0], ncbi, debug)
        if (rank3 == "superkingdom"):
            return (tax_name3.lower())

######################################################################
def name2taxid(taxids, ncbi):
    """Converting taxon names to taxids
    Created by Joe R. J. Healey; Nick Youngblut
    Original code.
    """
    
    # If names were provided in taxid list, convert to taxids
    #taxids = taxids.replace('"', '').replace("'", '').split(',')
    
    new_taxids = []
    for taxid in taxids:
        try:
            new_taxids.append(ncbi.get_name_translator([taxid])[taxid][0])
        except KeyError:
            try:
                new_taxids.append(int(taxid))
            except ValueError:
                print ()
                print (colored("** ERROR: Cannot convert to taxid: " + taxid, 'red'))
                exit()

    return new_taxids[0]

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
    db = sys.argv[1]
    taxid = sys.argv[2]
    option = sys.argv[3]
    debug=True

    ## print
    if debug:
        debug_message('arguments\n', "yellow")
        debug_message('db:' + db, "yellow")
        if taxid.isdigit():
            debug_message('taxid: ' + str(taxid) , "yellow")
        else:
            debug_message('taxid: ' + taxid, "yellow")
        debug_message('option:' + option, "yellow")

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

    string_info = [str(int) for int in info]
    #HCGB.functions.main_functions.printList2file("./out_file.txt", string_info)

    if option == 'unravel':
        ## assuming all belong to same superkingdom
        group_obtained = get_superKingdom(string_info[0], ncbi, debug)
        
        dict_entries = get_GenBank_ids(db, string_info, 10, debug, group_given=group_obtained)
    
############################################################
if __name__== "__main__":
    main()

