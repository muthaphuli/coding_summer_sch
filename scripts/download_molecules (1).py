# -*- coding: utf-8 -*-
"""
Dr Krishna K Govender
University of Johannesburg
06 February 2024
Script to run through a list of compound ids and download the SDF entries 
as well as write out the compound properties to a csv file using pubchempy 
"""

import pubchempy as pcp
import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

# Read the CID from file
df = pd.read_csv('reduced_list.txt', names=['cid'])

# Molecule properties
comp = []
formula = []
weight = []
smiles = []
iupac = []

# Check if sdfs directory exists
# If not create it
if os.path.exists('sdfs'):
    print('sdfs directory already exists')
else:
    os.mkdir('sdfs')

for index, rows in df.iterrows():
    comp += pcp.get_compounds(str(rows['cid']),'cid')
    print('Downloading:',str(rows['cid']))
    pcp.download('SDF', 'sdfs/'+str(rows['cid'])+'.sdf', str(rows['cid']), 'cid', overwrite=True)
    
print('Download completed')

for compound in comp:
    formula.append(compound.molecular_formula)
    weight.append(compound.molecular_weight)
    smiles.append(compound.isomeric_smiles)
    iupac.append(compound.iupac_name)
    
df['formula'] = formula
df['weight'] = weight
df['smiles'] = smiles
df['iupac'] = iupac

df.to_csv('pubchem_data.csv', sep=',', index=False)

# Check if imagess directory exists
# If not create it
if os.path.exists('images'):
    print('images directory already exists')
else:
    os.mkdir('images')

# Check if 3D directory exists
# If not create it
if os.path.exists('3D'):
    print('3D directory already exists')
else:
    os.mkdir('3D')
    
for index, rows in df.iterrows():
    m = Chem.MolFromSmiles(str(rows['smiles']))
    img = Draw.MolToImage(m)
    img.save('images/'+str(rows['cid'])+'.png')
    m2 = Chem.AddHs(m)
    AllChem.EmbedMolecule(m2)
    AllChem.MMFFOptimizeMolecule(m2)
    print(Chem.MolToMolBlock(m2),file=open('3D/'+str(rows['cid'])+'.mol','w+'))