#!/usr/bin/env python3
# VERSION HISTORY
# v: 0.0.1
# Written by: Hans Mueller Paul, June 2018
#                           NOTES


# IMPORT LIBRARIES
import sys
import argparse
import gffutils
from os import path
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

# DEFINING VARIABLES
parser = argparse.ArgumentParser()

parser.add_argument('-f', required=True, dest='f', 
                    help='path to input file in FASTA format'
                    )
parser.add_argument('-g', required=True, dest='f', 
                    help='path to input file in GFF format'
                    )
# parser.add_argument('-n', '--number', metavar='', dest='n', type=int, default=5,
#                    help='number of candidate primer pairs to pick, default = 5'
#                    )
parser.add_argument('-o', '--output', dest='o', default='output.fasta',
                    help='path to output file'
                    )
parser.add_argument('-v', '--verbose', action='store_true',
                    help='prints each step of each iteration (for debugging)'
                    )

args = parser.parse_args()




## extract all gene/transcript(mRNA) sequences
def parent_sequence(type):
    for parent in db.features_of_type(type):
        p_seq = parent.sequence(myFASTA)
        p_seq = Seq(p_seq, generic_dna)
        if parent.strand == '-':
            p_seq = p_seq.reverse_complement()
        print('>' + parent.id + '\n' + p_seq)

## extract all cds/exon/intron sequences under the same transcript id
def child_seq(type):
    for t in db.features_of_type('mRNA', order_by='start'): # or mRNA depending on the gff
        print('>' + t.id + '_' + type)
        # print(t.sequence(myFasta))
        seq_combined = ''
        for i in db.children(t, featuretype=type, order_by='start'): # or exon/intron
            seq = i.sequence(myFASTA, use_strand=False)  # use_strand doesn't work in 0.8; have to revcomp
            seq_combined += seq
        seq_combined = Seq(seq_combined, generic_dna)
        if t.strand == '-':
            seq_combined = seq_combined.reverse_complement()
        for i in range(0, len(seq_combined), 60):
            print(seq_combined[i:i+60])


def main():
    if args.verbose:
        print(args)

    # IMPORTING FASTA FILE
    with open('./Files/args.g', 'r') as f:
        myGFF = f.read()
    with open('./Files/args.f', 'r') as f:
        myFASTA = f.read()
    db = gffutils.create_db(myGFF, ':memory:', merge_strategy="create_unique", keep_order=True)


    # WRITE TO OUTPUT FILE
    with open('temp.fasta', 'w') as f2:
        f2.write(f"""

        """)

# CONFIRMATION MESSAGE
    print(f'The output file has been generated at {args.o}')

if __name__ == '__main__':
    main()