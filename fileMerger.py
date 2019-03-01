#!/usr/bin/env python3
"""
Takes all files in a given directory that follow a certain pattern and merges those together in to one big file

CFM_pipeline Part 7:
Last part of the pipeline combines all the output files to one big file again. Uses the output of spectraDataEditor.py
"""
import re;
import os;
import sys;

"""
Main function and runs everything in one go.
filePath = path to the files you wish to merge
fileName = the pattern of the file names "ether_131_part_" for example.
"""
def main(filePath, fileName):
	foundFiles = [f for f in os.listdir(filePath) if re.search(re.escape(fileName) + r'[0-9]{2}_output_MS2LDA_ready.mgf',f)];
	print(foundFiles)
	namePattern = foundFiles[0].split("/")[-1].split(".")[0];
	print(namePattern)
	newName = re.sub(r'[0-9]{1,5}_part_[0-9]{2}', 'full', namePattern);
	newFile = filePath + newName + ".mgf";
	outputFile = open(newFile,"w");
	for aFile in foundFiles:
		for line in open(filePath + aFile,"r"):
			outputFile.write(line);
	outputFile.close();
		
if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2]);
