#!/usr/bin/env python3
import sys;
from rdkit import Chem
from rdkit.Chem import AllChem
"""
This script takes a file containing ID smile combinations and neutralises the smiles (if possible) and then writes the
ID, neutralised smiles and old smiles to a new file

made by: Rutger Ozinga
30/10/2018 
"""
#----------------------------------------------------------------------------------------
""" code made by Hans de Winter, obtained from www.rdkit.org/docs-beta/Cookbook.html 
Comments added by Rutger Ozinga"""

"""creates a tuple of tuples containing the Neutralisation Reactions"""
def _InitialiseNeutralisationReactions():
	patts= (
		# Imidazoles
		('[n+;H]','n'),
		# Amines
		('[N+;!H0]','N'),
		# Carboxylic acids and alcohols
		('[$([O-]);!$([O-][#7])]','O'),
		# Thiols
		('[S-;X1]','S'),
		# Sulfonamides
		('[$([N-;X2]S(=O)=O)]','N'),
		# Enamines
		('[$([N-;X2][C,N]=C)]','N'),
		# Tetrazoles
		('[n-]','[nH]'),
		# Sulfoxides
		('[$([S-]=O)]','S'),
		# Amides
		('[$([N-]C=O)]','N'),
		)
		#Return a list of the molecules created from the Smarts string and Smiles String for each entry in patts
	return [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y,False)) for x,y in patts]

_reactions=None
"""takes a smiles string and reactions which is empty for some reason and uses the build in
 Neutralisation patterns to neutralize the given smile string"""
def NeutraliseCharges(smiles, reactions=None):
	global _reactions
	if reactions is None:
		if _reactions is None:
			_reactions=_InitialiseNeutralisationReactions()
	reactions=_reactions
	mol = Chem.MolFromSmiles(smiles)
	replaced = False
	for i,(reactant, product) in enumerate(reactions):
		while mol.HasSubstructMatch(reactant):
			replaced = True
			rms = AllChem.ReplaceSubstructs(mol, reactant, product)
			mol = rms[0]
	#returns the changes Smiles and true if the smile was neutralized or the original smile and false if it was not neutralized.
	if replaced:
		return (Chem.MolToSmiles(mol,True), True)
	else:
		return (smiles, False)
"""End of code made by Hans de Winter"""
#---------------------------------------------------------------------------------------

"""send the given file containing smiles through the Neutralisation method. adds the output to lineList"""
def smilesNeutralize(filePath):
	lineList = [];
	for line in open(filePath,"r"):
		splitLine = line.split(",");
		smile = splitLine[1].strip();
		(molSmiles, neutralised) = NeutraliseCharges(smile)
		#if molecule was neutralized
		if neutralised:
			lineList.append(splitLine[0] + "," + molSmiles + "," + smile + "\n");
		#if molecule not neutralized
		else:
			lineList.append(splitLine[0] + "," + smile + "\n");
	return lineList;

"""writes the content of lineList to a new file. This file will contain the ID, neutralized Smiles, and original Smiles string."""
def makeNeutralizedSmileFile(lineList, finalOutput):
	outputFile = open(finalOutput,"w");
	for line in lineList:
		outputFile.write(line);
	outputFile.close();
    
"""
Main Method calls all other methods
requires the path to a the folder that contains the csv file with ID smile combinations 
"""
def main(filePath):
	#path to the file that needs to be ran through cfm_id
	smileFile = filePath;
	#path to the new File. Split on / remove the last entry (old file Name) and then put it back together.
	splitPath = filePath.split(".");
	outPath = splitPath[0] + "_neutralized.csv"
	lineList = smilesNeutralize(filePath);
	makeNeutralizedSmileFile(lineList,outPath);
	return outPath;
    
if __name__ == "__main__":
	main(sys.argv[1]);

