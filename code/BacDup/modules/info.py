#!/usr/bin/env python3
##############################################################
## Jose F. Sanchez & Alba Moya                              ##
## Copyright (C) 2020-2021                                  ##
##############################################################
'''
Created in March 2021
@author: Jose F. Sanchez-Herrero


This module contains help information for multiple options

'''
from termcolor import colored

import HCGB.sampleParser
import BacDup

##########################
def run_info(options):
    
    ## project help
    if (options.help_project):
        project_help()
        #sampleParser.help_project()
        exit()

    ## help_format option
    if (options.help_format):
        sampleParser.help_format()
        exit()

    ## help_input option
    if (options.help_input):
        BacDup.modules.input_parser.help_input()
        exit()

    
    print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))

##########################
def project_help():
    print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))
    
