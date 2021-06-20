#!/usr/bin/env python3
##############################################################
## Jose F. Sanchez & Alba Moya                              ##
## Copyright (C) 2020-2021                                  ##
##############################################################
from Bio.Cluster import Tree
from pickle import TRUE
'''
Created on 3 dic. 2020
@author: alba

Modified in April 2021
@author: Jose F. Sanchez-Herrero
'''
import os
import sys
import pandas as pd
import numpy as np
import scipy
from termcolor import colored
from Bio import SeqIO, Seq
from collections import defaultdict 

import argparse
from argparse import ArgumentParser

import HCGB
from HCGB.functions.aesthetics_functions import debug_message
import HCGB.functions.time_functions as HCGB_time

import BacDup
import BacDup.scripts.functions as BacDup_functions

## set words to exclude regarding phages
words2exclude_phages = ['integrase', 'terminase',
                     'phage',  'lysin', 'endolysin', 'holin', 
                     'capsid', 'tail', 'bacteriophage', 
                     'portal', 
                     'tapemeasure', 'tape measure', 
                     'baseplate', 'base plate',
                     'virion', 'antirepressor', 'excisionase',
                     'Cro-like represessor', 'CI-like repressor'
                     'rIIA lysis', 'rI lysis', 'rIIB lysis',
                     'head decoration', 'HNH endonuclease', 'single-stranded DNA-binding protein']


################################################################################
def get_dup_stats(sample, dup_annot_df, annot_table, debug):
    '''Generate some statistics for the duplicated analysis
    '''
    
    data2add = pd.DataFrame(columns=BacDup_functions.columns_dup_table())
    
    ################################################
    ## add info for each grouping of duplicates
    ################################################
    total_dupGroups = len(dup_annot_df["dup_id"].unique())
    total_dups = dup_annot_df.shape[0]
    total_prots = annot_table.shape[0]
    total_dupGroups_pseudoFree = len(dup_annot_df["dup_id_pseudo_free"].unique())
    total_dupGroups_pseudoMobileFree = len(dup_annot_df["dup_id_mobile_free"].unique())

    ################################################
    ## get some information from duplicated genes
    ################################################
    # eg. transposases
    transposases_count = dup_annot_df[ dup_annot_df['product'].str.contains('transposase')].shape[0]
    
    # phage associated
    dup_annot_df_tmp = dup_annot_df.copy()
    for w in words2exclude_phages:
        dup_annot_df_tmp = dup_annot_df_tmp[ ~(dup_annot_df_tmp['product'].str.contains(w, regex=False ))]

    phage_count = dup_annot_df_tmp.shape[0]
    
    # hypothetical
    hypothetical_count = dup_annot_df[ dup_annot_df['product'].str.contains('hypothetical')].shape[0]
    
    # pseudo
    pseudo_count = dup_annot_df[ dup_annot_df['pseudo']==True].shape[0]
    
    # etc
    
    ## fill dataframe
    data2add.loc[sample, "total_prots"] = total_prots
    data2add.loc[sample, "n_groups_all"] = total_dupGroups
    data2add.loc[sample, "n_dups_all"] = total_dups
    data2add.loc[sample, "n_groups_pseudoFree"] = total_dupGroups_pseudoFree
    data2add.loc[sample, "n_groups_pseudoMobFree"] = total_dupGroups_pseudoMobileFree
    data2add.loc[sample, "transpo_count"] = transposases_count
    data2add.loc[sample, "phage_count"] = phage_count
    data2add.loc[sample, "hypothetical_coun"] = hypothetical_count
    data2add.loc[sample, "pseudo"] = pseudo_count
   
    print(f"""
    \t\t+ Sample: {sample} 
    \t\t {total_prots} total proteins annotated

    \t\t {total_dups} proteins duplicated 
    \t\t Duplicate proteins annotation (non-exclusive):
    \t\t  |- {phage_count} phage-associated proteins
    \t\t  |- {transposases_count} transposases proteins
    \t\t  |- {hypothetical_count} hypothetical proteins
    \t\t  |- {pseudo_count} pseudogenes

    \t\t Number of groups of duplicates including: """)
    
    ################################################
    ## create statistics for each group 
    ################################################
    dupGroups_pseudoFree = dup_annot_df[dup_annot_df["dup_id_pseudo_free"]>0]
    dupGroups_pseudoMobileFree = dup_annot_df[dup_annot_df["dup_id_mobile_free"]>0]
    
    groups_duplicates = {"All duplicates              " : ["all", dup_annot_df],
                         "Pseudogene Free             " : ["pseudoFree", dupGroups_pseudoFree],
                         "Pseudo, transpo & phage Free" : ["pseudoMobFree", dupGroups_pseudoMobileFree] }

    for name, df_dup in groups_duplicates.items():
        list_entries = df_dup[1].groupby(["dup_id"]).count()['count_dups'].to_list()
        list_entries = [x for x in list_entries if x==x] ## remove NaNs
        list_entries = [int(x) for x in list_entries] ## convert to numbers
        
        ## get distribution statistics
        ## get median, SD duplicates/group
        ## get biggest group
        d = scipy.stats.describe(list_entries)
        df = pd.DataFrame([d], columns=d._fields)
        #"nobs", "minmax", "mean", "variance", "skewness" , "kurtosis",
        # d[0]    d[1]       d[2]    d[3]        d[4]        d[5]
        
        ## convert to strings to save
        list_entries = [str(x) for x in list_entries]    
        mean = "{0:.3g}".format(d[2])
        counts_dups = len(df_dup[1]["dup_id"])
        n_tot = d[0]
        n_min = d[1][0]
        n_max = d[1][1]
        
        ## fill dataframe
        data2add.loc[sample, "n_groups_" + df_dup[0]] = n_tot
        data2add.loc[sample, "n_dups_" + df_dup[0]] = counts_dups
        data2add.loc[sample, "n_min_" + df_dup[0]] = n_min
        data2add.loc[sample, "n_max_" + df_dup[0]] = n_max
        data2add.loc[sample, "n_mean_" + df_dup[0]] = mean
        data2add.loc[sample, "list_" + df_dup[0]] = ":".join(list_entries)
        
        ## print in screen
        print(f"""\t\t  |- {name} = {n_tot} groups (n={counts_dups}; max={n_max}; min={n_min}; mean={mean})""")

    print()
    return(data2add)

################################################################################
def get_dupannot(sample, blast_results_df, annot_table_file, debug):
    '''Get an annotation information for duplicated proteins'''
    
    ## debug messages
    if debug:
        debug_message('get_dupannot function', 'yellow')
        debug_message('blast_results_df: ', 'yellow')
        print (blast_results_df)
        debug_message('annot_table_file: ' + annot_table_file, 'yellow')
        
    #get duplicated protein list
    qseqid = list(blast_results_df["qseqid"])
    sseqid =list(blast_results_df["sseqid"])
    
    qseqid.extend(sseqid) ## Fix me
    prot_id = list(set(qseqid))
    
    ## gets annotation
    annot_table = HCGB.functions.main_functions.get_data(annot_table_file, ',', "index_col=0")
    
    #get filtered_annot table
    filtered_annot = annot_table.loc[prot_id]

    ###################################
    ## Get duplicate groups    
    ###################################
    # 1st round
    relations_dict = defaultdict(list) 
    for index, row in blast_results_df.iterrows():
        relations_dict[row['qseqid']].append(str(row['sseqid']))
    
    ## debug messages
    if debug:
        debug_message('relations_dict: ', 'yellow')
        print (relations_dict)
        
    ## 2nd round
    new_relations_dict = defaultdict(list)
    dups=0
    for key, value in relations_dict.items():
        stop=False
        for dup_id, new_value in new_relations_dict.items():
            if key in new_value:
                stop=True
        if not stop:
            for key2, value2 in relations_dict.items():
                if (key == key2):
                    continue
                else:
                    if (key2 in value): 
                        for i in value2: 
                            if i not in value: 
                                value.append(str(i))
            dups += 1
            value.append(str(key))
            new_relations_dict[str(dups)] = value
            
            ## debug messages
            if debug:
                print()
                debug_message("dups: " + str(dups), 'yellow')
                debug_message("key: " + key, 'yellow')
                print(value)
                print()
                
    ## debug messages
    if debug:
        debug_message('new_relations_dict: ', 'yellow')
        print (new_relations_dict)

    ## Create data
    filtered_annot['dup_id'] = ""
    #df_data = pd.DataFrame(columns=('index', 'dup_id'))
    for dup_id, new_value in new_relations_dict.items():
        for i in new_value:
            filtered_annot.loc[i, 'tmp_dup_id'] = dup_id

    ## debug messages
    if debug:
        debug_message('filtered_annot: ', 'yellow')
        HCGB.functions.main_functions.print_all_pandaDF(filtered_annot)
        
    #filtered_annot = filtered_annot.drop(columns='index')

    ## Check group of duplicates
    ## group & count
    filtered_annot["count_dups"] = filtered_annot.groupby("tmp_dup_id")["tmp_dup_id"].transform("count")
    
    ## Some proteins might be orphans if for some reason they do not fulfill cutoffs with
    ## all members of a duplicated group. That might create some groups with only 1 protein.
    filtered_annot = filtered_annot[filtered_annot["count_dups"]>1]
    
    ## reset dup_ids
    df_grouped = filtered_annot.groupby('tmp_dup_id')
    dup_id_count=0
    dup_annot_df_fixed = pd.DataFrame()
    for group, df_group in df_grouped:
        dup_id_count += 1
        df_group['dup_id'] = dup_id_count
        dup_annot_df_fixed = dup_annot_df_fixed.append(df_group)

    ## add a column reflecting number of members with no PSEUDOGENES
    
    ## remove old dup_id column
    del dup_annot_df_fixed['tmp_dup_id']
    
    ######################################
    ## Filter pseudogenes
    ######################################
    print ("+ Pseudogenes would be used in the analysis and flagged appropriately...")
    
    ### debug
    #HCGB.functions.main_functions.print_all_pandaDF(dup_annot_df_fixed)
    pseudo_free = dup_annot_df_fixed.copy()
    pseudo_free = pseudo_free[pseudo_free['pseudo'] != True]

    ## debug messages
    if debug:
        debug_message('Dimension dataframe: ', 'yellow')
        print("Before:" + str(dup_annot_df_fixed.shape))
        print("After:" + str(pseudo_free.shape))

    ## Check group of duplicates without pseudo
    ## group & count
    pseudo_free["count_dups_pseudo_free"] = pseudo_free.groupby("dup_id")["dup_id"].transform("count")
    
    ## Some proteins might be orphans if for some reason they do not fulfill cutoffs with
    ## all members of a duplicated group. That might create some groups with only 1 protein.
    pseudo_free = pseudo_free[pseudo_free["count_dups_pseudo_free"]>1]
    
    ## reset dup_ids for pseudo_free
    dup_id_count2=0
    pseudo_free['dup_id_pseudo_free'] = ""
    df_grouped2 = pseudo_free.groupby('dup_id')
    pseudo_free2 = pd.DataFrame()
    for group, df_group2 in df_grouped2:
        dup_id_count2 += 1
        df_group2['dup_id_pseudo_free'] = dup_id_count2
        pseudo_free2 = pseudo_free2.append(df_group2)

    ## debug messages
    if debug:
        debug_message('pseudo_free annotation: ', 'yellow')
        print(pseudo_free2)
    
    ######################################
    ## Filter phages and transposases
    ######################################
    print ("+ Genes annotated as Mobile elements (phages and transposases) would be used in the analysis and flagged appropriately...")
    
    ## copy elements
    mobile_free = pseudo_free2.copy()

    # debug messages
    if debug:
        debug_message('Dimension dataframe: ', 'yellow')
        print("Before:" + str(dup_annot_df_fixed.shape))
        print("After:" + str(mobile_free.shape))

    ## get words to exclude
    words2exclude = words2exclude_phages + ['transposase']
    for w in words2exclude:
        before_count = mobile_free.shape[0]
        mobile_free = mobile_free[ ~(mobile_free['product'].str.contains(w, regex=False ))]
    
        ## TODO: save statistic into dataframe?
        
        # debug messages
        if debug:
            print("--------------------------")
            debug_message('Remove: ' + w, 'yellow')
            after_count = mobile_free.shape[0]
            print ("%s elements left " %after_count)
            print ("%s elements removed " %(before_count-after_count))
            
    
    ## Check group of duplicates without pseudo
    ## group & count
    mobile_free["count_dups_mobile_free"] = mobile_free.groupby("dup_id")["dup_id"].transform("count")
    
    ## Some proteins might be orphans if for some reason they do not fulfill cutoffs with
    ## all members of a duplicated group. That might create some groups with only 1 protein.
    mobile_free = mobile_free[mobile_free["count_dups_mobile_free"]>1]
    
    ## reset dup_ids for pseudo_free
    dup_id_count3=0
    mobile_free['dup_id_mobile_free'] = ""
    df_grouped3 = mobile_free.groupby('dup_id')
    mobile_free2 = pd.DataFrame()
    for group, df_group3 in df_grouped3:
        dup_id_count3 += 1
        df_group3['dup_id_mobile_free'] = dup_id_count3
        mobile_free2 = mobile_free2.append(df_group3)

    ## debug messages
    if debug:
        debug_message('mobile_free annotation: ', 'yellow')
        print(mobile_free2)

    ###################################
    ## Merge all duplicates groups
    ###################################
    
    ## debug messages
    if debug:
        debug_message('dup_annot_df_fixed: ', 'yellow')
        print (dup_annot_df_fixed)

    ## merge columns
    frames2merge = [dup_annot_df_fixed, 
                    pseudo_free2[['dup_id_pseudo_free', 'count_dups_pseudo_free']],
                    mobile_free2[['dup_id_mobile_free', 'count_dups_mobile_free']]
                    ]
    df_merge_all = pd.concat(frames2merge, axis=1, join='outer')
    
    ## get some statistics
    data2add = get_dup_stats(sample, df_merge_all, annot_table, debug)

    ## return         
    return(df_merge_all, data2add)

################################################################################
def create_blast_results(sample, fasta_file, outdir, debug):
    '''Creates BLAST results for each fasta vs. itself'''
    
    #phr is the header file, pin is the index file, psq is the sequence file
    
    ## debug messages
    if debug:
        debug_message('create_blast_results function call:', 'yellow')
        debug_message('sample: ' + sample, 'yellow')
        debug_message('fasta_file: ' + fasta_file, 'yellow')
        debug_message('outdir: ' + outdir, 'yellow')
    
    ## output file
    raw_blast = os.path.abspath(os.path.join(outdir, "BLAST_raw_results.tsv"))

    ## timestamps 
    db_timestamp = os.path.join(outdir, '.db_success')
    search_timestamp = os.path.join(outdir, '.blast_success')
        
    if (not HCGB.functions.files_functions.is_non_zero_file(search_timestamp)):

        ## get binaries
        (makeblastdb_exe, blastp_exe) = BacDup.modules.config.get_exe('BLAST', debug)
        makeblastdb_exe = "/usr/bin/makeblastdb" 
        blastp_exe = "/usr/bin/blastp"
        
        ## check if db is indexed already
        db_path_name = os.path.join(os.path.abspath(outdir), sample + '_db')
        if (not HCGB.functions.files_functions.is_non_zero_file(db_timestamp)):
            ## generate blastdb for genome
            HCGB.functions.blast_functions.makeblastdb(db_path_name, fasta_file, makeblastdb_exe, 'prot') # HCGB function    
        
            ## print time stamp
            HCGB_time.print_time_stamp(db_timestamp)
        
        else:
            print (colored("\t+ BLAST database already available for sample %s [%s]" %(sample, read_time), 'green'))
            
        ## create blastp outfile
        HCGB.functions.blast_functions.blastp(blastp_exe, raw_blast, db_path_name, fasta_file, 1) # HCGB function

        ## print time stamp
        HCGB_time.print_time_stamp(search_timestamp)
    else:
        read_time = HCGB_time.read_time_stamp(search_timestamp)
        print (colored("\t+ Duplicate search already available for sample %s [%s]" %(sample, read_time), 'green'))
            
    return (raw_blast)

################################################################################
def check_annot_table(annot_table, file, format, debug):
    '''Check annotation table provided matches the BLAST or protein fasta file provided'''
    
    ## read annotation
    annotation_table = pd.read_csv(annot_table, sep=",", index_col=0 )
    #annotation_table = annotation_table.drop("Unnamed: 0",axis=1)
    
    ## debug messages
    if debug:
        debug_message('check_annot_table function:', 'yellow')
        debug_message("list(annotation_table.columns)", 'yellow')
        print(list(annotation_table.columns))
        debug_message("BacDup_functions.columns_annot_table()", 'yellow')
        print(BacDup_functions.columns_annot_table())
    
    ## skip if not OK
    if not (list(annotation_table.columns) == BacDup_functions.columns_annot_table()):
        print('ERROR: Annotation table does NOT match desired input format')
        return(False)
    
    ## read BLAST sequence results
    if format=="blast":
        ## TODO
        print()
    
    ## protein fasta file
    elif(format=="fasta"):
        with open (file) as in_handle:
            ref_recs = SeqIO.to_dict(SeqIO.parse(in_handle, "fasta"))

        protein_IDs = list(ref_recs.keys()) ## get sequence headers
        index_annot = list(annotation_table.index) ## get annotation labels
        
        ## debug messages
        if debug:
            debug_message("list(annotation_table.index).sort()", 'yellow')
            print(index_annot)
            debug_message("list(ref_recs.keys()).sort()", 'yellow')
            print(protein_IDs)
        
        if (index_annot == protein_IDs):
            return(True)
        else:
            return(False)
        
################################################################################
def filter_BLAST(raw_blast, pident, evalue, bitscore, percent):
    '''Filter BLAST tabular results
    Created on 3 dic. 2020
    @author: alba
    '''    
    
    #fill aln_pct columns
    raw_blast["aln_pct_qlen"] = (raw_blast["length"]/raw_blast["qlen"]*100).round(2)
    raw_blast["aln_pct_slen"] = (raw_blast["length"]/raw_blast["slen"]*100).round(2)

    #filter mirror pairs
    del_mirror =pd.DataFrame(np.sort(raw_blast[["qseqid", "sseqid"]], axis=1))
    raw_blast = raw_blast[~del_mirror.duplicated()]
          
    ## filter using different parameters
    filter_pident= raw_blast["pident"] >= pident
    filter_evalue = raw_blast["evalue"] <= evalue
    filter_bitscore = raw_blast["bitscore"] >= bitscore
    filter_alnpct_slen = raw_blast["aln_pct_slen"] >= percent
    filter_alnpct_qlen = raw_blast["aln_pct_qlen"] >= percent
    filter_id = raw_blast["qseqid"] != raw_blast["sseqid"]

    filtered_results = raw_blast[filter_pident & filter_evalue & filter_bitscore & filter_alnpct_qlen & filter_alnpct_slen & filter_id]
    
    return(filtered_results)

################################################################################
def get_data(sample, file_data, format, out_folder, debug):
    '''Function to get BLAST results'''
    
    file_data = os.path.abspath(file_data)
        
    ## debug messages
    if (debug):
        debug_message('dup_searcher.get_data:', 'yellow')
        debug_message('file_data:' + file_data, 'yellow')
        debug_message('format:' + format, 'yellow')
    
    ## check file is readable
    BacDup_functions.file_readable_check(file_data)
    
    ## parse accordingly
    if (format=='blast_raw'):
        ## FIXME
        raw_blast = pd.read_csv(file_data, sep="\t", header = None, names=BacDup_functions.columns_rawBLAST_table())
                
    elif (format=='fasta'):
        
        raw_blast = create_blast_results(sample, file_data, out_folder, debug)
        raw_blast = pd.read_csv(raw_blast, sep="\t", header = None, names=BacDup_functions.columns_rawBLAST_table())
        
    return (raw_blast)

################################################################################
def filter_data(sample, file_data, format, pident, evalue, percent, bitscore, out_folder, debug):
    
    '''
    Get duplicated proteins and generate a filtered and sort dataframe depending on arguments values
    if a BLAST results text file is provided, it should have next columns names.    
    '''
    ## get BLAST results and filter
    raw_blast = get_data(sample, file_data, format, out_folder, debug)
    filtered_results = filter_BLAST(raw_blast, pident, evalue, percent, bitscore)

    #sort by aln_pct_qlen (desc), evalue(asc), bitscore(desc)
    by_alnpct = filtered_results.sort_values(["aln_pct_qlen", "evalue", "bitscore"],
                                       ascending=[False, True, False])
    return(by_alnpct)

###############################################################
def help_options():
    print ("\nUSAGE: python %s file format pident--evalue--percent--bitscore Debug[True/False]\n"  %os.path.realpath(__file__))

############################################################
def main():
    ## this code runs when call as a single script

    ## control if options provided or help
    if len(sys.argv) > 1:
        print ("")
    else:
        help_options()
        exit()
        
    file_data = os.path.abspath(sys.argv[1])
    format = sys.argv[2]
    filters = sys.argv[3].split('--')
    debug = sys.argv[4]
        
    filter_data(file_data, format, pident, evalue, percent, bitscore, '.', debug)

############################################################
if __name__== "__main__":
    main()
