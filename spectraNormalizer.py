#!/usr/bin/env python3

"""
Takes the mgf spectra files made with cfm Id and normalizes them on a scale from 0 to 900. Writes the output to new file.
requires the path to the spectra and a name pattern (for example, ether_131_part_). 

CFM_pipeline Part 5:
Fourth part of the pipeline is dependant on the output of CFMrunner.py. Requires the spectra output files. next part in the pipeline is 
tandemMS_Merger.py.
Made by: Rutger Ozinga 
Last edit: 10/10/2018
"""

import math;
import re;
import os;

"""
out of a file path it saves all the spectra in to a list of lists. 
each list contains all the lines of a single spectra.
Also runs the normalization so in one go all the data is normalized and put in to a List of lists.
This is then returned.
"""
def makeSpectraList(filePath):
    spectraList = [];
    spectra = [];
    for line in open(filePath):
        #end of a spectra is indicated with END IONS so when it appears 
        #all data stored in the the list is added to the spectraList
        #and the list is emptied for the next spectra
        if line.startswith("END IONS"):
            if spectra != []:
                spectra.append(line.strip());
                weightList, normList = normalizeSpectra(spectra);
                adv_normList = advancedNormalization(normList);
                #zips the list of weights and normalizations together.
                spectraZip = list(zip(weightList,adv_normList));
                #sends the spectra (list of all the lines in the spectra) and the spectraZip to the method remakeSpectra.
                spectraList.append(remakeSpectra(spectra,spectraZip));
                spectra = [];
        else: 
            spectra.append(line.strip());
    return spectraList;

"""
applies normal normalisation on the entered list containing the data of one spectra.
does it by calculating the mean and standard deviation, then it subtracts the mean each of the values.
then it divides the mean subtracted by the std. 
"""    
def normalizeSpectra(spectraList):
    total = 0;
    intensityList = [];
    weightList = []
    # for line in the peaks of the spectra list first 4 lines are headers and the last line isn't a peak either
    for line in spectraList[4:-1]:
        intensity = float(line.split()[1]);
        weightList.append(line.split()[0]);
        intensityList.append(intensity);
        total += intensity;
        n = len(spectraList[4:-1]);
    mean = total/n;
    #for x in intensityList subtract the mean for each value 
    meanSubbedIntensities = [x - mean for x in intensityList];
    #for each value in meanSubbedIntensities square it and sum them all up. 
    #then take the square root of that resulting value and save it as the std.
    std = math.sqrt(sum([x**2 for x in meanSubbedIntensities])/n);
    # in case the lenght of the spectra peak list is 1 (so only 1 peak) it will just return the mean.
    if(n == 1):
        normIntensities = [mean];

    else:
        #for x in meanSubbedIntensities divided it by the std
        normIntensities = [x/std for x in meanSubbedIntensities];
    return weightList, normIntensities;
"""
takes the normalized values and further normalizes it by putting it on a 0-900 scale. first it sets the values of 
everything above 2 and below -2 to 2 or -2 respectivly.
then it multiplies the values by 225. making the values go from -450 to 450
then by adding 450 the min value is now 0 and the highest values 900.
expects a list of normalizedIntensities.
"""
def advancedNormalization(normIntensities):
    for x in range(len(normIntensities)):
        if normIntensities[x] > 2:
            normIntensities[x] = 2.0;
        elif normIntensities[x] < -2:
            normIntensities[x] = -2.0;
    magicInducedList = [x*225 for x in normIntensities];
    moreMagic = [x + 450 for x in magicInducedList];
    return moreMagic;
"""
takes the original spectra and the normalized one and creates a new spectra taking 
the headers from the original and the peaks from the normalized list of peaks.
"""
def remakeSpectra(originalSpectra, normalizedSpectra):
    newSpectra = [];
    for i in range(len(originalSpectra)):
        if(i < 3):
            newSpectra.append(originalSpectra[i]);
        elif(3 < i < len(originalSpectra)-1): 
            newSpectra.append(str(normalizedSpectra[i-4][0]) + " " + str(normalizedSpectra[i-4][1]));
        else:
            newSpectra.append(originalSpectra[i]);
    return newSpectra;
    
"""
takes a list of lists containing the normalized spectra and puts them in to a new files 
requires a list of lists and a file path.
"""
def writeNewFile(newFile,normSpecList):
    normSpecFile = open(newFile,'w');
    for spectra in normSpecList:
        for line in spectra:
            normSpecFile.write(str(line) + "\n");
    normSpecFile.close();

"""
runs the main method.
requires the path to a the result folder that contains the mass spectra 
and the common part of the name of the file "ethers_131_part_" for example.
"""
def main(resultPath,fileName):
    filePath = resultPath;
    #for file in directory if file ends with 2 numbers followed by _output.mgf save the file in a list.
    #match example   "smiles_1000_part_01_output.mgf"
    foundFiles = [f for f in os.listdir(filePath) if re.search(re.escape(fileName) + r'[0-9]{2}_output.mgf',f)];
    for aFile in foundFiles:
        splitName = aFile.split(".");
        #reformat the old name with normalized in it 
        #output example given "smiles_1000_part_01_output.mgf" --> "smiles_1000_part_01_output_normalized.mgf"
        newName = splitName[0] + "_normalized." + splitName[1];
        spectraList = makeSpectraList(filePath + "" + aFile);
        newFilePath = filePath + "" + newName;
        writeNewFile(newFilePath,spectraList);

if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2]);

