#!/usr/bin/env python
"""
This script validates structures, converts to canonical representation

Usage example:
```
python rdkit_validator.py canonize input.txt output.txt
```

Second example provided below expects both files having canonical SMILES 
representations:
```
python rdkit_validator.py intersect reference.txt generated.txt
```
It doesn't converts anything to SMILES, simply compares these files and writes 
to stdout line with number of intersection and total smiles strings found in 
second file.
"""

import sys, argparse
from rdkit import Chem

def process_file(filename):
    f = filename
    for line in f:
        molecule = Chem.MolFromSmiles(line.strip())
        if not molecule is None:
            yield Chem.MolToSmiles(molecule)


def canonization_func(args):
    g = args.output_filename
    for m in process_file(args.filename):
        g.write("%s\n" % m)


def count_intersections_func(args):
    # for now use memory-consuming representation with sets
    a = set([l.strip() for l in args.filename if l.strip() != '']) 
    b = set([l.strip() for l in args.filename2 if l.strip() != ''])
    print("Intersects: %s, total: %s"% (len(a & b), len(b)))


    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    canonization_parser = subparsers.add_parser("canonize", 
        help="canonize SMILES list and save to file")
    canonization_parser.add_argument("filename", 
        help="file with list of formulae in SMILES, one per line", 
        type=argparse.FileType('r'))
    canonization_parser.add_argument("output_filename", 
        help="output file - to save list of canonical SMILES, one per line", 
        type=argparse.FileType('w'))
    canonization_parser.set_defaults(func=canonization_func)
    
    intersections_parser = subparsers.add_parser("intersect", 
        help="count number of SMILES in second file which present in 1st file")
    intersections_parser.add_argument("filename",
        help = "reference list of canonical smiles structures", 
        type=argparse.FileType('r'))
    intersections_parser.add_argument("filename2", 
        help="filename with smiles to check against reference", 
        type=argparse.FileType('r'))
    intersections_parser.set_defaults(func=count_intersections_func)
    
    args = parser.parse_args()
    args.func(args)
    