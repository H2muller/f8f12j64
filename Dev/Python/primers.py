#!/usr/bin/env python3
# VERSION HISTORY
# v: 0.0.1
# Written by: Hans Mueller Paul, March 2018
#                           First iteration of the code with set values for all variables
#                           code is functional until determination of TM and GC content
#                           list filtering for GC and TM not working
# v: 0.0.2
# Edited by: Hans Mueller Paul, May 2018
#                           Added Argument Parser section to the code and associated modifications
#                           fixed issue with GC / TM filtering. Primer pairing based on TM not working
#                           output file generation code not written
# v: 0.0.3
# Edited by: Hans Mueller Paul, June 2018
#                           [notes on version]

# IMPORT LIBRARIES
# import platform
# import sys
import datetime
import itertools
import argparse
import math
import statistics
from tabulate import tabulate
# import os
from os import path

# DEFINING VARIABLES
parser = argparse.ArgumentParser()
parser.add_argument('-i', required=True, metavar='', dest='i', help='path to input file in FASTA format')
parser.add_argument('-n', '--number', metavar='', dest='n', type=int, default=5, help='number of candidate primer pairs to pick, default = 5')
parser.add_argument('-e', '--extension', metavar='', dest='e', type=int, default=100, help='number of bases from the start and end of the sequence to look for primers in, default =100')
parser.add_argument('-s', '--short', metavar='', dest='s', type=int, default=20, help='shortest acceptable primer, default = 20')
parser.add_argument('-l', '--long', metavar='', dest='l', type=int, default=30, help='longest acceptable primer, default = 30')
parser.add_argument('-m', '--mintemp', metavar='', dest='m', type=float, default=55, help='min Tm in celsius, default = 55')
parser.add_argument('-x', '--maxtemp', metavar='', dest='x', type=float, default=62, help='max Tm in celsius, default = 62')
parser.add_argument('-M', '--mingc', metavar='', dest='M', type=float, default=40, help='min GC%, default = 40')
parser.add_argument('-X', '--maxgc', metavar='', dest='X', type=float, default=60, help='max GC%, default = 60')
parser.add_argument('-D', '--tmdiff', metavar='', dest='D', type=float, default=0.5, help='accepted TM difference to form primer pair')
parser.add_argument('-o', '--output', metavar='', dest='o', default='output.fasta', help='path to output file in FASTA format')
parser.add_argument('-v', '--verbose', action='store_true', help='prints each step of each iteration (for debugging)')

args = parser.parse_args()

# UNUSED VARIAVLES (copied from the perl software by dr. Hudson, for reference)
# -e <exponent-mantissa pair> (min e-value for primer BLAST check, default = 1)
# -b <boolean> (if TRUE (-b T specified), BLAST step is skipped)
# -t <integer> (if integer specified, tile primer pairs with gap specified across sequence, default = one primer at each end)
# -o <integer> (overlap for tiled fragments, ignored if -t not specified, default = 100)
# -p <integer> max successive homopolymer bases tolerated in primer, default = 4
# -q <boolean> (if TRUE specified, primers must contain all upper-case bases)
# -a <integer> position within string to start looking for primers (default = 0)

def main():
    if args.verbose:
        print(args)

    # IMPORTING FASTA FILE
    if args.verbose:
        print('Item exists: ' + str(path.exists(args.i)))
    if path.exists(args.i):
        with open(args.i, 'r') as f:
            contents = f.read()
            seqlist = contents.split()
        if args.verbose:
            print('seqlist is {}'.format(type(seqlist)))
        nametag = seqlist[0]
        seq = ''.join(seqlist[1:])
        if args.verbose:
            print(f'file name is {nametag} and sequence is {seq}')
    else: print('file does not exist')

    # DEFINING REVERSE COMPLIMENT
    alt_map = {'ins':'0'}
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    def revcomp(input):
        for k,v in alt_map.items():
            compseq = seq.replace(k,v)
        bases = list(compseq) 
        bases = reversed([complement.get(base,base) for base in bases])
        bases = ''.join(bases)
        for k,v in alt_map.items():
            bases = bases.replace(v,k)
        return bases
    
    # DEFINING PRIMER LISTS
    def get_primer(input_string):
        length = len(input_string[:args.e])
        prim_list_1 = [input_string[i:j + 1] for i in range(length) for j in range(i+args.s,i+args.l)]
        if args.verbose:
            print('prim_list_1 is a {}'.format(type(prim_list_1)))
        return prim_list_1
    if args.verbose:
        print(get_primer(seq))
    def fwd(input):
        return get_primer(input)
    def rvs(input):
        seq_rc = revcomp(input)
        return get_primer(seq_rc)
    fwd_primers = fwd(seq)
    rvs_primers = rvs(seq)
    if args.verbose:
        print(f'List of forward primers: {fwd_primers}')
        print(f'List of reverse primers: {rvs_primers}')

    # DEFINING Primer GC CONTENT
    def GC(input):
        gc_count = input.count('g') + input.count('c') + input.count('G') + input.count('C')
        gc_fraction = float(gc_count)/len(input)
        return 100 * gc_fraction

    # DEFINING PRIMER MELTING TEMPERATURES
    def TM(input):
        N = len(input)
        melt = 81.5 + 0.41*(GC(input)) - (675/N)
        return melt

    # DELETING FORWARD PRIMERS THAT DO NOT FALL WITHIN %GC AND TM LIMITS

    if args.verbose:
        print(f'Number of forward primers prior to GC verification: {len(fwd_primers)}')
    # fwd_primers = list(filter(GCrange,fwd_primers))
    for sequence in fwd_primers:
        if GC(sequence) < args.M: 
            fwd_primers.remove(sequence)
        elif GC(sequence) > args.X: 
            fwd_primers.remove(sequence)
    if args.verbose:
        print(f'Number of forward primers after to GC verification: {len(fwd_primers)}')
        print(f'Number of forward primers prior to TM verification: {len(fwd_primers)}')
    # fwd_primers = list(filter(TMrange,fwd_primers))
    for sequence in fwd_primers:
        if TM(sequence) < args.m: 
            fwd_primers.remove(sequence)
        elif TM(sequence) > args.x: 
            fwd_primers.remove(sequence)
    if args.verbose:
        print(f'Number of forward primers after to TM verification: {len(fwd_primers)}')

    # DELETING REVERSE PRIMERS THAT DO NOT FALL WITHIN %GC AND TM LIMITS
    if args.verbose:
        print(f'Number of reverse primers prior to GC verification: {len(rvs_primers)}')
    for sequence in rvs_primers:
        if GC(sequence) < args.M: 
            rvs_primers.remove(sequence)
        elif GC(sequence) > args.X: 
            rvs_primers.remove(sequence)
    if args.verbose:
        print(f'Number of reverse primers after to GC verification: {len(rvs_primers)}')
        print(f'Number of reverse primers prior to TM verification: {len(rvs_primers)}')
    for sequence in rvs_primers:    
        if TM(sequence) < args.m    : 
            rvs_primers.remove(sequence)
        elif TM(sequence) > args    .x: 
            rvs_primers.remove(sequence)
    if args.verbose:
        print(f'Number of reverse primers after to TM verification: {len(rvs_primers)}')

    # FORMING PRIMER PAIRS
    if args.verbose:
        print(f"Pairing primers based on TM (primers in a pair will have a maximum TM difference of {args.D} degrees)")
    
    PrimerPairs = list(itertools.product(fwd_primers,rvs_primers))
    if args.verbose:
        print(f"Total number of primer pairs: {len(PrimerPairs)}")
    Pairs = [i for i in PrimerPairs if math.isclose(TM(i[0]),TM(i[1]),abs_tol=args.D)]
    if args.verbose:
        print(f"Number of primer pairs with matching TM: {len(Pairs)}")
    Temps = [(TM(i[0]),TM(i[1])) for i in Pairs]
    if args.verbose:
        print(f"TMs for each primer pair: {Temps}")
    MatchingTemp = list(zip(Pairs, Temps))
    if args.verbose:
        print(f"matched pairs: {MatchingTemp}")
    primertable = tabulate(MatchingTemp, headers=['Primer pair','TM of primer pair'])

    # WRITE TO OUTPUT FILE
    with open(args.o, 'w') as f2:
        f2.write(f"""
        Target gene: {nametag}

        Genetic sequence:
        {seq}

        List of primers with TMs:
        {primertable}
        
        Generated with $software_name in {datetime.datetime.now()}
        Please cite $citation_quote
        """)

# CONFIRMATION MESSAGE
    print(f'The output file has been generated at {args.o}')

if __name__ == '__main__':
    main()