'''
Created on 5 ene. 2021

@author: alba
'''
parser = ArgumentParser(prog='dupAnnotation',
                        formatter_class=argparse.RawDescriptionHelpFormatter,
                        description="Get an annotation file with duplicated protein on genome.")

annotparser = parser.add_argument_group('Annotation parser options named arguments')
annotparser.add_argument("-a", "--annot_file", metavar="", help="Annotation file: genbank or GFF.")
annotparser.add_argument("-r", "--ref_file", metavar="", help="Genome references FASTA file.")
###
blastoptions = parser.add_argument_group('BLAST options named arguments')
blastoptions = parser.add_argument_group('BLAST options named arguments')
blastoptions.add_argument("-b", "--blast_folder", metavar="", help="BLAST binary folder. **Note blast_folder=/usr/bin is set by default**")    
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

arg = c(parser.parse_args(), 
arg_dict = vars(arg)