#!/usr/bin/env python3

"""
Takes the the mgf output of cfm_predict thats normalized and merged and adds the ID, InchiKey, smile and IUPAC fields to each spectra.

CFM_pipeline Part 6:
Takes the output of steps 2 and 5 and puts them together in new files that end with MS2LDA_Ready.mgf these files can already be put in MS2LDA.
However 1 more step will merge the splitted files back together. (INSERT NAME OF NEXT STEP)
Made by: Rutger Ozinga 
Last edit: 10/10/2018
"""

import os;
import sys;
import re;

"""
Turns the given spectra file in to a a list of lists. each list contains 1 spectra (from Begin IONS to End IONS)
first variable is th path to the spectra file (Mgf format).
"""
def makeSpectraList(spectraPath):
    spectraList = [];
    spectra = [];
    for line in open(spectraPath):
        #end of a spectra is indicated with END IONS so when it appears 
        #all data stored in the the list is added to the spectraList
        #and the list is emptied for the next spectra
        if line.startswith("END IONS"):
            if spectra != []:
                spectra.append(line);
                #sends the spectra (list of all the lines in the spectra) and the spectraZip to the method remakeSpectra.
                spectraList.append(spectra);
                spectra = [];
        else: 
            spectra.append(line.strip());
    return spectraList;
    
"""
Adds the info from the spectra infolist (id,smile and inchikey) to the spectra in the spectra list.
each spectra finds the title which in the case of cfm id is in location 4 for the spectra (so 3 position in a list).
Finds the spectra id in the given spectrainfo dict, and merges the info. Also adds a IUPAC field when the gets added to the database 
that way it can be added to the MS2LDA output. In case the smiles string was charged 2 additional fields are added another smile and inchiKey
field for the neutral smiles string and its corresponding inchiKey.
spectraList = list of lists with each list containing a spectra.
spectraInfo = a dictionary containing the NP-DB ids as key and a list of smile, inchikey combinations as values
"""
def addSpectraInfo(spectraList,spectraInfo):
    editedSpectraList = [];
    for spectra in spectraList:
        editedSpectra = [];
        cutUpTITLE = spectra[3].split(";");
        print(cutUpTITLE);
        structure_id = cutUpTITLE[0].split("=")[1];
        print(structure_id);
        ID = "ID=" + structure_id;
        newTitle = "TITLE=" + ";".join(cutUpTITLE[1:]);
        spectraData = spectraInfo[structure_id];
        if len(spectraData) == 2:
            smile = "SMILES=" + spectraData[0];
            inchiKey = spectraData[1];
            for x in range(len(spectra)):
                if x < 3:
                    if x == 1:
                        editedSpectra.append("IUPAC={}".format("Not_Added"));
                        editedSpectra.append(ID);
                        editedSpectra.append(newTitle);
                    editedSpectra.append(spectra[x]);
                elif x == 3:
                    editedSpectra.append(smile);
                    editedSpectra.append(inchiKey);
                else:
                    editedSpectra.append(spectra[x]);
            editedSpectraList.append(editedSpectra);
        else:
            #the spectraData consists of neutralSmile, smile, neutralInchiKey, inchiKey
            neutralSmiles = "NEUTRAL_SMILES=" + spectraData[0];
            smile = "SMILES=" + spectraData[1];
            neutralInchiKey = "NEUTRAL_" + spectraData[2];
            inchiKey = spectraData[3];
            #adds the new info between the already existing info in the spectra.
            for x in range(len(spectra)):
                if x < 3:
                    if x == 1:
                        #IUPAC. Add the IUPAC Key here in the format to add it to the MGF output file.
                        editedSpectra.append("IUPAC={}".format("Not_Added"));
                        editedSpectra.append(ID);
                        editedSpectra.append(newTitle);
                    editedSpectra.append(spectra[x]);
                elif x == 3:
                    editedSpectra.append(smile);
                    editedSpectra.append(neutralSmiles);
                    editedSpectra.append(inchiKey);
                    editedSpectra.append(neutralInchiKey);
                else:
                    editedSpectra.append(spectra[x]);
            editedSpectraList.append(editedSpectra);
    return editedSpectraList;
            
"""
Turns the content of the given path in to a dictionary with the first entry of the line as key and the second and third as a list values
first vaiable contains the path to the file. file should contain lines with 3 tab separated values in this case ID, smiles string, inchikey
"""
def spectraInfoDict(infoPath):
    spectraInfo = {};
    for line in open(infoPath):
        splittedLine = line.split();
        spectraInfo[splittedLine[0]] = splittedLine[1:];
    return spectraInfo;

"""
writes the info from the first variable in to the given file path.
first variable is a list of lists 
second variable is a path
"""
def writeFile(editedSpectra, newFile):
    editSpecFile = open(newFile,'w');
    for spectra in editedSpectra:
        for line in spectra:
            editSpecFile.write(str(line) + "\n");
    editSpecFile.close();

"""
Main function runs the rest of the script.
smilePath = path to the smile data
resultPath = path to the spectra result files
fileName = beginning part of the file names (ethers_131_part_) 
"""
def main(smilePath,resultPath,fileName):
    nameSplit = fileName.split("_");
    datafile = "_".join(nameSplit[:-3]);
    infoPath = smilePath + datafile + "_dataFile.txt";
    spectraPath = resultPath;
    #for file in directory if file ends with 2 numbers followed by _output.mgf save the file in a list.
    #match example   "smiles_1000_part_01_output.mgf"
    foundSpectraFiles = [f for f in os.listdir(spectraPath) if re.search(re.escape(fileName) + r'[0-9]{2}_output_normalized_merged.mgf',f)];
    foundSpectraFiles.sort();
    spectraInfo = spectraInfoDict(infoPath);
    for x in range(len(foundSpectraFiles)):
        spectraFilePath = spectraPath + foundSpectraFiles[x];
        spectraList = makeSpectraList(spectraFilePath);
        editedSpectra = addSpectraInfo(spectraList,spectraInfo);
        newFile = spectraFilePath.replace("normalized_merged","MS2LDA_ready");
        os.system("rm {}".format(spectraPath + foundSpectraFiles[x]));
        writeFile(editedSpectra,newFile);
    
if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2],sys.argv[3]);
