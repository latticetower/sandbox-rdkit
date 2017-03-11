#!/usr/bin/env python
"""
This scripts represents pipeline which do the following:
1. reads data from sdf file and converts it to fingerprint,
and saves it to table with solubility information
2. The results are processed with xgboost and trained model is saved somewhere.
3. given the results from p.1 and p.2, trained model tries to predict solubility
for structures in test dataset

All parts are independent and might be implemented as a separate scripts.
"""
import argparse, os
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

sdf_file_name = "solubility_2007.sdf"


def get_fingerprint(molecule, sep=","):
    """
    This is moved to separate function to change in 1 place when necessary
    returns string (probably ready for saving to .csv)
    warning!!! don't use GetMorganFingerprint(...)
    """
    return sep.join([str(x) for x in AllChem.GetMorganFingerprintAsBitVect(molecule, 2)])
    
    
def load_database(filename, sep=","):
    if not os.path.exists(filename):
        print("sdf database couldn't be found at '%s'" % filename)
        exit(1)
    suppl = Chem.SDMolSupplier(filename)
    for mol in suppl:
        yield sep.join([mol.GetProp('EXPT'), get_fingerprint(mol, sep=sep)])
        
    
if __name__=="__main__":
    with open('logS_data.csv', 'w') as f:
        for x in load_database(sdf_file_name):
            f.write("%s\n" % x)
