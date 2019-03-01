#!/usr/bin/env python3
import sys;
import re;
import os;

"""send the given file containing smiles through the Neutralisation method. adds the output to lineList"""
def smilesSplit(filePath):
	lineList = [];
	for line in open(filePath,"r"):
		splitLine = line.split(" ");
		originalSmile = splitLine[1].strip();
		splittedOriginalSmile = originalSmile.split(".");
		longestSmilePart = max(splittedOriginalSmile, key=len);
		lineList.append(splitLine[0] + " " + longestSmilePart + "\n");
	return lineList;

"""writes the content of lineList to a new file. This file will contain the ID, neutralized Smiles, and original Smiles string."""
def makeEditedSmileFile(lineList, finalOutput):
	outputFile = open(finalOutput,"w");
	for line in lineList:
		outputFile.write(line);
	outputFile.close();
   
"""
get all the files in a given file path that corresponds to a given pattern.
smilePath = path to the smile file
filePattern = patern that matches to the first part of the file
"""
def getFileList(smilePath,filePattern):
	pattern = re.escape(filePattern) + r'[0-9]{2}.txt';
	fileList = [f for f in os.listdir(smilePath) if re.search(pattern,f)];
	#sort it so the files go from 00 to 09;
	fileList.sort();
	return fileList;
    
"""
Main Method calls all other methods
requires the path to a the folder that contains the csv file with ID smile combinations 
"""
def main(filePath, fileName):
	#path to the new File. Split on / remove the last entry (old file Name) and then put it back together.
	fileList = getFileList(filePath,fileName);
	#list of running jobs
	for x in range(len(fileList)):
		newFilePath = filePath + fileName;
		lineList = smilesSplit(filePath + fileList[x]);
		makeEditedSmileFile(lineList,fileList[x]);

    
if __name__ == "__main__":
	main(sys.argv[1],sys.argv[2]);
