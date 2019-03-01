#!/usr/bin/env python

"""
Merges the output of cfm-id (mgf format), which contains 3 different energy level, in to one spectra.

CFM_pipeline Part 6:
takes the output of spectraNormalizer and merges the normalized spectra of each molecule (cfm creates 3 spectra for each molecule on 3 energy levels.)
It creates new files containing the normalized merged spectra. next step in the pipeline is spectraDataEditor.py.

Made by: Rutger Ozinga 
Last edit: 10/10/2018
"""

import os;
import re;


energyList = [10,20,40];

"""
Merges a dict containing weights as key and a list of intensities as value. in to a weight with a single intensity.
also combines weights that are similar enough (currently put on a range of 1 of eachother) and removes the entries with a intensity below the user given cutoff value.
teruns a list of sorted weights (low to high).
"""
def spectraMerge(spectraDict,cutoff):
    mergedList = [];
    #iterate over all keys in the dict. In this case the weights
    keyList = spectraDict.keys();
    keyList = sorted(keyList);
    prevWeight = 0
    for weight in spectraDict.keys():
        specSum = 0
        #loop through the intensities in the list.
        for intensity in spectraDict[weight]:
            specSum += intensity;
        #divide the intensity by the length of the entry to normalize the values.
        mergedSpec = [weight, round(specSum/len(spectraDict[weight]))]; 
        #filter out low intensities with a user given cutoff value
        if mergedSpec[1] > cutoff:
            mergedList.append(mergedSpec);
        specSum = 0;
    mergedList.sort();
    return mergedList;
    
    
    
"""
this method takes as input a list lists containing the data of 3 spectra
each list in the list contains the data of one spectra in mgf format.
It then combines these 3 spectra in to 1 spectra. First it takes the first 4 lines of the first spectra.
This contains the header of the spectra, the only difference between these 4 lines in the 3 spectra is the line 
"Energyxx" this tells what energy level was used in cfm-id. since we combine it, its changed to EnergyCombined along with the 3 energy levels used.
after these four lines the next lines till it reaches the end will contain a peak with its weight and intensity.
These are added to a dictionary. Key is the weight and the value is a list containing the intensities. 
All peaks with the same weight are combined and the intensities added to the list corresponding to the weight

it returns the dict.

"""
def combineSpectraSet(spectraSet,cutoff):
    combinedSpectra = [];
    spectraDict = {};
    for x in range(3):
        if x == 0:
            #for the first spectra it will add the first 4 lines of the mgf file to the new combined spectra.
            for headerLine in range(4):
                #change the ID line from Energy0 to EnergyCombi
                if headerLine == 3:
                    newEnergy = 'EnergyCombined {0}eV {1}eV {2}eV'.format(energyList[0],energyList[1],energyList[2]);
                    replaced = re.sub('Energy\d', newEnergy, spectraSet[x][headerLine]);
                    combinedSpectra.append(replaced);
                else:
                    combinedSpectra.append(spectraSet[x][headerLine]);
                    
        for y in range(4,len(spectraSet[x])-1):
            splitSpectra = spectraSet[x][y].split();
            if splitSpectra[0] in spectraDict:
                # append the new weight to the existing array in the dict
                spectraDict[float(splitSpectra[0])].append(float(splitSpectra[1]));
            else:
                # create a new array of intensities in the dict with the given weight as key
                spectraDict[float(splitSpectra[0])] = [float(splitSpectra[1])];
    #calls the method spectraMerge which will merge the dict in to a List of lists of weight spectra combinations.
    mergedList = spectraMerge(spectraDict, cutoff);
    #loops through the list of lists and adds the spectra to a list containing strings.
    for spectra in mergedList:
        combinedSpectra.append(str(spectra[0]) + " " + str(spectra[1]));
    combinedSpectra.append("END IONS");
    return combinedSpectra;

"""
takes a list of spectra and each step it takes 3 spectra (belonging to low medium and high voltages)
then combines them in a list and sends the list of 3 spectra to the function combineSpectraSet.
returns the combined colllection of the spectra (so 3 spectra combined in to 1 done for all the spectra in spectaList).
"""      
def createSpectraSets(spectraList,cutoff):
    combinedCollection = [];
    listLen = len(spectraList);
    spectraSet = [];
    #from 0 to X with steps of 3  example : 0,3,6,9,etc
    for i in range(0,listLen,3):
        spectraSet.append(spectraList[i]);
        spectraSet.append(spectraList[i+1]);
        spectraSet.append(spectraList[i+2]);
        #add the output of combineSpectra to list of combined spectra.
        combinedCollection.append(combineSpectraSet(spectraSet,cutoff));
        #empty spectraSet for the next 3 spectra in the loop.
        spectraSet = [];
    return combinedCollection;
    
"""
reads to a given path adds line to list ill it reaches END IONS then it 
adds the list to a spectra List and continues with the next.
requires a filepath to run.
"""
def createSpectraList(filePath):
    spectraList = [];
    spectra = [];
    for line in open(filePath):
        if line.startswith("END IONS"):
            if spectra != []:
                spectra.append(line.strip());
                spectraList.append(spectra)
                spectra = [];
        else: 
            spectra.append(line.strip());
    return spectraList;
    
"""
takes a list of lists containing the combined spectra and puts them in to a new files 
requires a list of lists and a file path.
"""
def writeNewFile(newFile,combiSpecCol):
    newSpecFile = open(newFile,'w');
    for spectra in combiSpecCol:
        for line in spectra:
            newSpecFile.write(line + "\n");
    newSpecFile.close();
    
"""
main method runs the other methods.
requires the path to a the result folder that contains the mass spectra 
and the common part of the name of the file "ethers_131_part_" for example.
""" 
def main(resultPath, fileName, optionDict):
    filePath = resultPath;
    foundFiles = [f for f in os.listdir(filePath) if re.search(re.escape(fileName) + r'[0-9]{2}_output_normalized.mgf',f)];
    for aFile in foundFiles:
        print("Merging: " + filePath + aFile);
        splitName = aFile.split(".");
        newName = splitName[0] + "_merged." + splitName[1];
        spectraList = createSpectraList(filePath + aFile);
        #contains all the combined spectra in a list of lists.
        combiSpecCol = createSpectraSets(spectraList,optionDict[cutoff_intensity]);
        newFilePath = filePath + "" + newName;
        writeNewFile(newFilePath,combiSpecCol);
        os.system("rm {}".format(filePath + aFile));
        
if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2]);

            
