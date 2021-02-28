'''
Created on 24 nov. 2020

@author: alba
'''

import os
import sys

def create_folder(path):
    """Create a folder directory 'path'. Returns path created."""

       
    #access_rights = 0o755 readable by all but write access by owner
    #access_rights = 0o777 (default) readable and writable by all
    # define the access rights absolutely necessarY?
    try:
        os.mkdir(path)
    except OSError:  
        print ("\tDirectory %s already exists\n" %path)
        return path
    else:  
        print ("Successfully created the directory %s\n" %path)
        

    return path

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print (__doc__)
        
       
        print ("## Usage create_folder: ")
        print ("python %s path\n" %sys.argv[0])
       

        sys.exit()

    create_folder(*sys.argv[1:])