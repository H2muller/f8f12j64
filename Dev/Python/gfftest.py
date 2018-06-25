#!/usr/bin/env python3

# Libraries
import gffutils
import os

# Function variables
with open('./Files/sbicolor.gff3', 'r') as f:
    pathtogff = f.read()
fn = gffutils.example_filename(pathtogff)
db = gffutils.create_db(fn, dbfn='test.db', keep_order=True)

def main():
    print(open(fn).read()) 


if __name__ == '__main__':
    main()