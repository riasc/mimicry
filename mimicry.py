from pyfaidx import Fasta
from pysam import FastaFile 
from icecream import ic
from Bio import SeqIO
from intervaltree import Interval, IntervalTree
import json
import sys
import configargparse
import random
import logging
from pathlib import Path


class Variants():

    def __init__(self, args):
        self.variants = {}

        # multiple fasta files
        self.args = args
        self.get_genome()

        # create output directory
#        make_dir(args.outdir)

        self.svnum = {'DEL': args.del_num, 
                 'INS': args.ins_num, 
                 'INV': args.inv_num, 
                 'DUP': args.dup_num}

        self.svlens = {'DEL': args.del_len, 
                  'INS': args.ins_len, 
                  'INV': args.inv_len, 
                  'DUP': args.dup_len}

        
        # simple variants
        self.generate_variants()
        self.change_genome()


#        print(self.variants)


    def get_genome(self):
        self.genome = Fasta(self.args.ref)
        # initiate interval tree for each chromosome 
        self.variants = {key: IntervalTree() for key in self.genome.keys()}
        logging.info('reference genome loaded: %s', self.args.ref)

    def change_genome(self):
        for chrom in self.variants.keys():
            for interval in self.variants[chrom]:
                variant = interval.data
                    
                

    # generates deletions, insertions, and inversions
    def generate_variants(self):
        for i in range(0, max(self.svnum.values())):
            if i < self.svnum['DEL']:
                self.place_variants('DEL', i)
            if i < self.svnum['INS']:
                self.place_variants('INS', i)
            if i < self.svnum['INV']:
                self.place_variants('INV', i)

    def place_variants(self, svtype, i):
        # check if the interval/variant is already present in the tree
        (chrom, start, end) = self.generate_random_interval(svtype)

        # check if the interval/variant is already present in the tree
        ovlp = self.variants[chrom].overlap(start, end)
        if len(ovlp) != 0: # overlapping, try other position
            for j in range(0, self.args.tries):
                (chrom, start, end) = self.generate_random_interval(svtype)
                ovlp = self.variants[chrom].overlap(start, end)
                if len(ovlp) == 0:
                    break

        # still overlapping
        if len(ovlp) != 0:
            mess = 'Could not find non-overlapping position'
            mess += 'for variant of type %s'
            logging.warning(mess, svtype)
            logging.warning('exceeded number of tries (%d)',
                            self.args.tries)
            logging.warning('skipping - %s variants placed', i)

        else:
            vcfentry = self.generate_vcf_entry(i, chrom, start, end, svtype)
            self.variants[chrom][start:end] = vcfentry


    def generate_random_interval(self, svtype):
        chrom = self.random_chromosome()
        start = random.randint(0, len(self.genome[chrom])-1)
        svlen = random.randint(int(self.svlens[svtype][0]),
                               int(self.svlens[svtype][1]))
        return (chrom, start, start + svlen - 1)


    def random_chromosome(self):
        return random.choice(list(self.genome.keys()))


    # def make_dir(path):
        # if not self.args.outdir.exists():
            # self.args.outdir.mkdir(parents=True, exist_ok=True)


    def generate_vcf_entry(self, i, chrom, start, end, svtype):
        vcfentry = {'CHROM': chrom, 'POS': str(start+1),
                    'ID': 'mimicry.'+svtype+'.'+str(i),
                    'REF': '.', 'ALT': '.', 'QUAL': '.',
                    'FILTER': 'PASS', 'INFO': '.', 'FORMAT': '.'}


        info = {'SVTYPE': svtype, 'END': str(end+1),
                'SVLEN': str(end-start+1)}

        if svtype == 'DEL':
            vcfentry['REF'] = self.genome[chrom][start:end].seq
            vcfentry['ALT'] = self.genome[chrom][start-1].seq
            info['SVLEN'] = str(-1 * (end-start+1))
        elif svtype == 'INS':
            vcfentry['REF'] = self.genome[chrom][start-1].seq
            vcfentry['ALT'] = self.genome[chrom][start:end].seq
        elif svtype == 'INV':
            vcfentry['REF'] = self.genome[chrom][start:end].seq
            vcfentry['ALT'] = self.genome[chrom][start:end].reverse.complement.seq


        vcfentry['INFO'] = ';'.join([key+'='+val for key, val in info.items()])


        return vcfentry




def main():
    parser = configargparse.ArgParser(default_config_files=['params.cfg'])
    parser.add('-c', '--config', required=True, 
               is_config_file=True, 
               help='configuration file for the simulation')

    parser.add('--ref', help='the reference genome to place the variants on',
               type=str)

    parser.add('--del_num', help='number of deletions to be simulated', 
               type=int)
    parser.add('--del_len', help='size of the deletions to be simulated',
               action='append')
    parser.add('--ins_num', help='number of deletions to be simulated', 
               type=int)
    parser.add('--ins_len', help='size of the deletions to be simulated',
               action='append')
    parser.add('--inv_num', help='number of deletions to be simulated', 
               type=int)
    parser.add('--inv_len', help='size of the deletions to be simulated',
               action='append')
    parser.add('--dup_num', help='number of deletions to be simulated', 
               type=int)
    parser.add('--dup_len', help='size of the deletions to be simulated',
               action='append')
    parser.add('--tries', help='number of tries to find a non-overlapping position',
               type=int)


    parser.add('--outdir', help='output directory to store the results',
               type=Path)

    args = parser.parse_args()


    # set up logging
    logging.basicConfig(filename=args.outdir / 'mimicry.log', 
                        level=logging.DEBUG)


    v = Variants(args)


main()



