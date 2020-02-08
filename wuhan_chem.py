#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 19:35:02 2020

@author: andyjiang
"""

from rdkit import Chem
from rdkit.Chem import AllChem

peptides = {
        
        'alanine' : 'N[C@@H](C)C(=O)O',
        'glycine' : 'NCC(=O)O',
        'aspartic acid' : 'N[C@@H](CC(=O)[O-])C(=O)O',
        'glutamic acid' : 'N[C@@H](CCC(=O)[O-])C(=O)O',
        'lysine' : 'N[C@@H](CCCC[NH3+])C(=O)O',
        'arginine' : 'N[C@@H](CCCNC(=[NH2+])N)C(=O)O',
        'histidine' : 'N[C@@H](CC1=CN=CN1)C(=O)O',
        'valine' : 'N[C@@H](C(C)C)C(=O)O',
        'isoleucine' : 'N[C@@H](C(CC)C)C(=O)O',
        'serine' : 'N[C@@H](CO)C(=O)O',
        'threonine' : 'N[C@@H](C(O)C)C(=O)O',
        'tyrosine' : 'N[C@@H](Cc1ccc(O)cc1)C(=O)O',
        'asparagine' : 'N[C@@H](CC(=O)N)C(=O)O',
        'glutamine' : 'N[C@@H](CCC(=O)N)C(=O)O',
        'tryptophan' : 'N[C@@H](CC1=CNc2ccccc21)C(=O)O',
        'phenylalanine' : 'N[C@@H](Cc1ccccc1)C(=O)O',
        'proline' : 'N1[C@@H](CCC1)C(=O)O',
        'methionine' : 'N[C@@H](CCSC)C(=O)O',
        'leucine' : 'N[C@@H](CC(C)C)C(=O)O',
        'cysteine' : 'N[C@@H](CS)C(=O)O'
        
        
        }


conversion = {
        'A' : 'alanine',
        'G' : 'glycine',
        'D' : 'aspartic acid',
        'E' : 'glutamic acid',
        'K' : 'lysine',
        'R' : 'arginine',
        'H' : 'histidine',
        'V' : 'valine',
        'I' : 'isoleucine',
        'S' : 'serine',
        'T' : 'threonine',
        'Y' : 'tyrosine',
        'N' : 'asparagine',
        'Q' : 'glutamine',
        'W' : 'tryptophan',
        'F' : 'phenylalanine',
        'P' : 'proline',
        'M' : 'methionine',
        'L' : 'leucine',
        'C' : 'cysteine'
            
        }

def letterToAA(letter):
    
    return conversion[letter]

#forms a polypeptide SMILES string from a Peptide Sequence and Ionizes it at pH 7.4
def polyPeptide(aaArr):
    
    #polypeptide String
    ppS = ""
    
    for i in range(len(aaArr) - 1):
        ppS = ppS + aaArr[i][0:len(aaArr[i])-1]
        
    ppS = ppS + aaArr[len(aaArr) - 1]
    
    
    
    if aaArr[0] == peptides['proline']:
        ppS = '[N1H2+]' + ppS[2:len(ppS)-1] + '[O-]'
        
    else:
        ppS = '[NH3+]' + ppS[1:len(ppS)-1] + '[O-]'
        
    
        
    return ppS

peptideFile = open("peptideFile.mol", "w")


"""

ala = 'N[C@@H](C)C(=O)O'
ala2 = '[NH3+][C@@H](C)C(=O)[O-]' 
test = '[' + ala[0:1] + 'H3+]' + ala[1:len(ala)-1] + '[' + ala[len(ala)-1:len(ala)] + '-]'
gly = 'N[C@@H]C(=O)O'

"""

"""

peptides = {
        
        'alanine' : 'N[C@@H](C)C(=O)O',
        'glycine' : 'NCC(=O)O',
        'aspartic acid' : 'N[C@@H](CC(=O)[O-])C(=O)O',
        'glutamic acid' : 'N[C@@H](CCC(=O)[O-])C(=O)O',
        'lysine' : 'N[C@@H](CCCC[NH3+])C(=O)O',
        'arginine' : 'N[C@@H](CCCNC(=[NH2+])N)C(=O)O',
        'histidine' : 'N[C@@H](CC1NC=NC1=)C(=O)O',
        'valine' : 'N[C@@H](C(C)C)C(=O)O',
        'isoleucine' : 'N[C@@H](C(CC)C)C(=O)O',
        'serine' : 'N[C@@H](CO)C(=O)O',
        'threonine' : 'N[C@@H](C(O)C)C(=O)O',
        'tyrosine' : 'N[C@@H](Cc1ccc(O)cc1)C(=O)O',
        'asparagine' : 'N[C@@H](CC(=O)N)C(=O)O',
        'glutamine' : 'N[C@@H](CCC(=O)N)C(=O)O',
        'tryptophan' : 'N[C@@H](CC1=CNc2ccccc21)C(=O)O',
        'phenylalanine' : 'N[C@@H](Cc1ccccc1)C(=O)O',
        'proline' : 'N1[C@@H](CCC1)C(=O)O',
        'methionine' : 'N[C@@H](CCSC)C(=O)O',
        'leucine' : 'N[C@@H](CC(C)C)C(=O)O',
        'cysteine' : 'N[C@@H](CS)C(=O)O'
        
        
        }


conversion = {
        'A' : 'alanine',
        'G' : 'glycine',
        'D' : 'aspartic acid',
        'E' : 'glutamic acid',
        'K' : 'lysine',
        'R' : 'arginine',
        'H' : 'histidine',
        'V' : 'valine',
        'I' : 'isoleucine',
        'S' : 'serine',
        'T' : 'threonine',
        'Y' : 'tyrosine',
        'N' : 'asparagine',
        'Q' : 'glutamine',
        'W' : 'tryptophan',
        'F' : 'phenylalanine',
        'P' : 'proline',
        'M' : 'methionine',
        'L' : 'leucine',
        'C' : 'cysteine'
            
        }

"""

#The polypeptide sequence of the wuhan coronavirus surface glycoprotein

wuhansg_file = open("wuhansg.txt", "r")
    
wuhan_sg = wuhansg_file.read().replace("\n","")

wuhan_sg_aaArr = []

for alpha in range(len(wuhan_sg)):
    lc = letterToAA(wuhan_sg[alpha:alpha+1])
    wuhan_sg_aaArr.append(peptides[lc])
    
wuhan_sg_Smiles = polyPeptide(wuhan_sg_aaArr)


mol = Chem.MolFromSmiles(wuhan_sg_Smiles)

mol = Chem.AddHs(mol)

#AllChem.EmbedMolecule(mol)
#AllChem.MMFFOptimizeMolecule(mol)

molBlock = Chem.MolToMolBlock(mol)

peptideFile.write(molBlock)

peptideFile.close()

"""

for aminoAcid, aaSmiles in peptides.items():
    try:
        testing = Chem.MolFromSmiles(aaSmiles)
        testing = Chem.AddHs(testing)
        
    except:
        print(aminoAcid, aaSmiles)
        
"""

"""

peptideFile = open("peptideFile.mol", "r")

pepLines = peptideFile.readlines()

print(pepLines[4][0:33].split())

peptideFile.close()

"""



"""

def letterToAA(letter):
    
    return conversion[letter]

#forms a polypeptide SMILES string from a Peptide Sequence and Ionizes it at pH 7.4
def polyPeptide(aaArr):
    
    #polypeptide String
    ppS = ""
    
    for i in range(len(aaArr) - 1):
        ppS = ppS + aaArr[i][0:len(aaArr[i])-1]
        
    ppS = ppS + aaArr[len(aaArr) - 1]
    
    if aaArr[0] == peptide['proline']:
        ppS = '[N1H2+]' + ppS[2:len(ppS)-1] + '[O-]'
        
    else:
        ppS = '[NH3+]' + ppS[1:len(ppS)-1] + '[O-]'
        
    return ppS
    
"""