#!/usr/bin/env python3
# Written by: Hans Müller Paul and Jacob Heldenbrand
#                           NOTES


# IMPORT LIBRARIES
import sys
import argparse
import gffutils
import re
from os import path
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

# DEFINING VARIABLES
parser = argparse.ArgumentParser()

parser.add_argument('-f', '--fasta', required=True, dest='f', 
                    help='path to input file in FASTA format'
                    )
parser.add_argument('-g', '--gff', required=True, dest='f', 
                    help='path to input file in GFF format'
                    )
parser.add_argument('-l', '--length', metavar='', dest='l', type=int, default=500,
                    help='length of flanking region for verification, default = 500'
                    )
parser.add_argument('-o', '--output', dest='o', default='output.fasta',
                    help='path to output file'
                    )
parser.add_argument('-v', '--verbose', action='store_true',
                    help='prints each step of each iteration (for debugging)'
                    )

args = parser.parse_args()


class Crispr:
    def __init__(self, sequence):
        self.sequence = sequence
        self.crispr = self.__search_for_crispr_target()
        self.sgrna = self.__design_sgrna_sequence()
        self.flanking = self.__flanking_region_for_sequencing()
    
    def __search_for_crispr_target(self):
        target = 'NGGN'
        crispr_target = [target.start() for target in re.finditer('?=%s' %(target),'%s' %(self))]
        return crispr_target
    
    def __design_sgrna_sequence(self):
        return

    def __flanking_region_for_sequencing(self):
        offset = args.l
        start = self-offset
        end = self+len(self.sgrna)+fpoffset
        # if start <=0:
        #     start == 0
        # if end >= len(sequence):
        #     end == len(sequence)
        genomic_region = sequence[start:end]
        return genomic_region


def get_sequence_by_type(type):
    for parent in featuredb.features_of_type(type):
        parent_seq = parent.sequence(myFASTA)
        parent_seq = Seq(parent_seq, generic_dna)
        if parent.strand == '-':
            parent_seq = parent_seq.reverse_complement()
        return ('>' + parent.id + '\n' + parent_seq)

def get_features_under_transcript(type):
    for transcript in featuredb.features_of_type('mRNA', order_by='start'): # or mRNA depending on the gff
        print('>' + transcript.id + '_' + type)
        seq_combined = ''
        for child in featuredb.children(transcript, featuretype=type, order_by='start'): # can also order by exon/intron
            child_seq = child.sequence(myFASTA, use_strand=False)  # use_strand is currently bugged, so I have it disabled
            seq_combined += child_seq
        seq_combined = Seq(seq_combined, generic_dna)
        if transcript.strand == '-':
            seq_combined = seq_combined.create_reverse_complement()
        for i in range(0, len(seq_combined), 60):
            print(seq_combined[i:i+60])



def main():
    if args.verbose:
        print(args)

    
    # IMPORT FASTA FILE AND CREATE DATABASE
    with open('./Files/args.g', 'r') as g:
        myGFF = g.read()
    with open('./Files/args.f', 'r') as f:
        myFASTA = f.read()
    featuredb = gffutils.create_db(myGFF, ':memory:', merge_strategy="create_unique",keep_order=True)
    if args.verbose:
        if len(featuredb) > 0:
            print('database successfully generated')

    # DETERMINE TARGET REGIONS (by position)


    # GENERATE TARGET SEQUENCES
    do search stuff
    generate sequence
    add to dictionary {gene_id:sequence}
    if bla.strand == '-':
        do revcomp, then stuff
    return dictionary


    # WRITE TO OUTPUT FILE
    with open('temp.fasta', 'w') as f2:
        f2.write(f"""

        """)

# CONFIRMATION MESSAGE
    print(f'The output file has been generated at {args.o}')

if __name__ == '__main__':
    main()