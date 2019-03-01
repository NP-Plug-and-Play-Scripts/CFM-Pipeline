#!/usr/bin/env python3

"""
This script turns a file of NP-ID and SMILES strings combinations in to a file that contains NP-ID, SMILES and InchiKey or NP-ID, neutral SMILES, SMILES, neutral InchiKey, Inchikey.
WARNING this requires molconvert a tool included in jchem. https://chemaxon.com/products/instant-jchem/download. 

CFM_pipeline step 3:
Second step in the pipeline. takes each of the splitted files and adds the inchi keys to it so they can later be added to the mgf files.
next script run is CFMrunner.py.

Made by: Rutger Ozinga 
Last edit: 22/11/2018
"""
import os;
import re;
import sys;
from rdkit import Chem;
import time

"""
Takes path and a file name (csv file) containing ids and smiles strings. 
I then with the help of rdkit turns the smiles in to their Kekule form. This prevents the error "Cannot process aromatic bonds" in molconvert.
filePath = path to the location of the file you wish to open.
"""
def kekuleSmilesMaker(smilesList):
	#new instance of SmilesFixes
	kekuleSmileList = [];
	count = 0
	for smiles in smilesList:
		count += 1;
		#run the fix_smiles method of smileFixer on the smile obtained frome the line.
		molecule = Chem.MolFromSmiles(smiles)
		Chem.Kekulize(molecule)
		#Chem.SanitizeMol(molecule,sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_KEKULIZE^Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
		fixedSmile = Chem.MolToSmiles(molecule, kekuleSmiles=True)
		kekuleSmileList.append(fixedSmile);
	return kekuleSmileList;

"""
creates 3 lists with the given file.
File should contain a ID, a neutral smiles string and an original smiles string.
"""
def createLists(path):
	idList = [];
	smileList = [];
	altSmileList = [];
	for line in open(path):
		line = line.strip();
		if path.endswith(".txt"):
			splittedLine = line.split(" ");
		elif path.endswith(".csv"):
			splittedLine = line.split(",");
		idList.append(splittedLine[0]);
		if len(splittedLine) > 2:
			smileList.append(splittedLine[1]);
			altSmileList.append(splittedLine[2]);

		else:
			smileList.append(splittedLine[1]);
			altSmileList.append(splittedLine[1]);
	return idList, smileList, altSmileList;

"""
creates a newFile from a list of smiles.
smileList = list of smiles strings
newPath = path for the new file to be placed
"""
def newFile(smileList,newPath):
	newFile = open(newPath,"w");
	for line in smileList:
		newFile.write(line + "\n");
	newFile.close();

"""
Runs the molconverter tool and takes the output file, reads it and puts the lines in a list and returns those.
The generated file gets deleted in the end.

"""
def createInchiKeys(molConvertPath,tempPath,smilePath):
	print(smilePath)
	#run commandline program molconvert to turn a file of smiles strings in to a file of inchikeys.
	os.system("{0} inchikey:SAbs {1} -o {2}".format(molConvertPath + "molconvert" ,smilePath,tempPath));
	inchiKeyList = [];
	for line in open(tempPath,"r"):
		inchiKeyList.append(line);
	os.system("rm {}".format(tempPath));
	return inchiKeyList;

"""
takes the id list (contains all the id's) smile list, alternative smile list, inchiKeyList, alternative inchikeyList and the path to 
the new file and writes the data in to the file
"""
def makeInchiSmileFile(idList,smileList,altSmileList,inchiKeyList, altInchiKeyList, finalOutput):
	finalOut = open(finalOutput,"w");
	print(len(idList), len(smileList),len(altSmileList), len(inchiKeyList), len(altInchiKeyList));
	for i in range(len(idList)):
		if smileList[i] != altSmileList[i]:
			#strip for the new inchi keys because molconvert realy likes to add line endings
			newLine = idList[i] + " " + smileList[i] + " " + altSmileList[i] + " " + inchiKeyList[i].strip() + " " + altInchiKeyList[i];
			finalOut.write(newLine);
		else:
			newLine = idList[i] + " " + smileList[i] + " " + inchiKeyList[i];
			finalOut.write(newLine);
	print("Data succesfully writen to a new file.")
	finalOut.close();

"""
main method runs the entire inchiKeyCreator. requires a path to molConvert, a path (path to the file directory) and a file name (name of the file)
file Number is optional depending on if you run it with the pipeline or not.  this makes it so the temporary files dont overwrite eachother.
"""
def main(molConvertPath,filePath, fileName, fileNumber = ""):
	start = time.time()
	print('working in ' + filePath + fileName);
	splitName = fileName.split(".");
	path = filePath + fileName;
	#new files to store the smiles in before they enter molconvert
	newSmilePath = filePath + "tempSmiles_" + fileNumber + ".txt";
	newNeutralSmilePath = filePath + "tempNeutralSmiles_" + fileNumber + ".txt";
	#location to store the inchikeys made with molconvert before adding them to the final file
	tempInchiPath = filePath + "tempInchiKeys_" + fileNumber + ".txt";
	tempNeutralInchiPath = filePath + "tempNeutralInchiKeys_" + fileNumber + ".txt";
	#new path for the output of the dataFile. Will contain the ID,SMILES,InchIKey, neutral Smiles and  neutral inchiKey.
	finalOutput = filePath + splitName[0] + "_dataFile.txt";
	#get a idList a smileList and altSmileList from the function createLists using the path to the smilesFile.
	idList, smileList, neutralSmileList = createLists(path);
	#turns the smileList and neutralSmileList in to kekule form.
	kekuleSmileList = kekuleSmilesMaker(smileList)
	neutralKekuleSmileList = kekuleSmilesMaker(neutralSmileList)
	#creates 2 new file swith the lists of kekuleSmiles.
	newFile(kekuleSmileList,newSmilePath);
	newFile(neutralKekuleSmileList, newNeutralSmilePath);
	#creates a list of inchiKeys with the createInchiKey method using the smiles.
	inchiKeyList = createInchiKeys(molConvertPath,tempInchiPath,newSmilePath);
    #creates a list of neutral inchiKeys with the createInchiKey method using the neutral smiles.
	neutralInchiKeyList = createInchiKeys(molConvertPath,tempNeutralInchiPath,newNeutralSmilePath);
	makeInchiSmileFile(idList, smileList, neutralSmileList, inchiKeyList, neutralInchiKeyList, finalOutput);
	os.system("rm {}".format(newSmilePath));
	os.system("rm {}".format(newNeutralSmilePath));
	end = time.time()
	print("Run took " + str(round((end - start) /60))+ " Minutes")
	
if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2],sys.argv[3]);

