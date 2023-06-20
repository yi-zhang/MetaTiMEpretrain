# function to process user input variables

import os
import sys
import argparse
import logging

#from metatimetrain import call_mec
#from metatimetrain import decompose
from metatimetrain import decompose_all

__doc__ = """"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('exprmatlist',help="list of scRNA expression files in h5ad or cell by gene txt")
    parser.add_argument('-f', '--module', nargs='+', default = ['decompose_all'], help = 'modules you want to run')
    
    parser.add_argument('-t', '--threads', default=4, type=int)
    parser.add_argument('-o','--output',default='./', type=str)

    args = parser.parse_args()

    #config output dir
    if not os.path.isdir(args.output): 
        os.makedirs(args.output)

    args.comptablefile = 'exprmatlist' #[putfile]

    # run the module requested
    if('decompose_all' in args.module):
        decompose_all.main(args.exprmatlist, args.threads)
    if('callmec' in args.module):
        # need : file of full component table.
        call_mec.main(args.comptablefile, args.output)
