from pyfaidx import Fasta
from pysam import FastaFile 
from icecream import ic
from Bio import SeqIO
from intervaltree import Interval, IntervalTree
import json
import sys
import configargparse

class Deletion():

    def __init__(self, args, no):
        for i in range(0,no):
            print(i)

class Variants():

    def __init__(self, args):
        self.variants = {}
        self.lens = {}

        self.genome = Fasta(args.ref)





def main():
    parser = configargparse.ArgParser(default_config_files=['params.cfg'])
    parser.add('-c', '--config', required=True, 
               is_config_file=True, 
               help='configuration file for the simulation')
    parser.add('--ref', help='the reference genome to place the variants on',
               type=str, action='append'))
    parser.add('--num_del', help='number of deletions to be simulated', 
               type=int)


    args = parser.parse_args()


    v = Variants(args)


main()



