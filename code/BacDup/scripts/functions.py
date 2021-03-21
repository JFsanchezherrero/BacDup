import HCGB

def pipeline_header(option):
    """
    Prints a common header for the pipeline including name, author, copyright and year.        
    """
    
    print ("\n")
    HCGB.functions.aesthetics_functions.print_sepLine("#", 70, False)
    
    print('#', '{: ^66}'.format("BacDup pipeline"), '#')
    print('#', '{: ^66}'.format("Jose F. Sanchez & Alba Moya Garces"), '#')
    print('#', '{: ^66}'.format("Copyright (C) 2020-2021"), '#')
    HCGB.functions.aesthetics_functions.print_sepLine("#", 70, False)
