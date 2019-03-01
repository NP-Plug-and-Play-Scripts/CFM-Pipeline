#!/usr/bin/env python3
"""
This script splits a given file in to smaller parts containing a equal amount of lines (or close to in case the total is uneven)
the cores is the amount of parts you want to split the file in to.

CFM_pipeline step 1:
This is the first script to be ran in the pipeline it takes the given .csv file and splits it in to equal parts. 
Done to reduce the time needed in the cfm_predict step. Next script in the pipeline is createInchiKeys.py.
Made by: Rutger Ozinga 
Last edit: 10/10/2018
"""

import sys;

#the amount of cores you wish to use.
cores = 10;
#general name of the files.
fileName = "";

"""
splits the given list in to a predivined number of parts and stores them in a dictionary
the way it divides the list is by taking the inner and outerr most enties 
of the length sorted List and adds them to the first entry.
It then moves up one spot at the start and down one spot a the end of the list and 
adds those to the next entry.
"""
def lineDivider(sortedList):
	global fileLength;
	fileLength = len(sortedList);
	fileNum = 1;
	partitionDict = {};
	for x in range(int(round(fileLength/2))):
		# if the fileNum is already in the dict append the existing list belonging to the key
		if fileNum in partitionDict.keys():
			#add the Xth entry of a list to the dict
			partitionDict[fileNum].append(sortedList[x]);
			#add the Xth entry at the end of a list to the dict 
			#example: list = [1,2,3,4,5] index 0 would be 1 from the start and index -1 would be 5 so to get that its -x-1 (-0 -1 = -1) 
			partitionDict[fileNum].append(sortedList[-x-1]);
		# if the fileNum is not a key in the dict add at with the values at the given index as a list of values.
		else:
			partitionDict[fileNum] = [sortedList[x],sortedList[-x-1]];
		#in case the number of cores has reached the max reset it to 1;
		if fileNum == cores:
			fileNum = 1;
		else:
			fileNum += 1;
	#in case the list has a uneven number of lines add the remaining line to the next file.
	if fileLength%2 != 0:
		partitionDict[fileNum].append(sortedList[int(fileLength/2)-1]);
	return partitionDict;
"""
takes  a file and adds all lines to a list. Then sorts the List based on length.
This method requires a fielPath.
"""
def fileLengthSorter(filePath):
	fileList = [];
	for line in open(filePath,"r"):
		fileList.append(line.replace(","," "));
	fileList = sorted(fileList,key=len);
	return fileList;

"""
given a dictionary containing numbers as keys and a list of lines it puts the content of the library
in individual files. So the values belonging to Key 1 will be added a file and the values of Key 2 to another file.
"""
def createFiles(partDict,smileFile,newPath):
	fileName = smileFile.split("/")[-1].split(".")[0];
	newFileName = fileName + "_" + str(int(round(fileLength/cores))) + "_part_";
	for part in partDict.keys():
		#formats the number in a way that instead of being 0,1,2,3 it will be 00,01,02,03
		partNum = "{0:02d}".format(part-1)
		filePart = open(newPath + newFileName + partNum + ".txt", "w"); 
		for line in partDict[part]:
			splitLine = line.split();
			filePart.write(splitLine[0] + " " + splitLine[1] + "\n");
	return newFileName;
"""
Main Method calls all other methods
requires the path to a the folder that contains the csv file with ID smile combinations 
"""
def main(filePath):
	print(filePath);
	#path to the file that needs to be ran through cfm_id
	smileFile = filePath;
	#path to the new File. Split on / remove the last entry (old file Name) and then put it back together.
	splitPath = filePath.split("/");
	newPath = "/".join(splitPath[:-1]) + "/";
	#length of the file. is calculated in lineDivider()
	fileLength = 0;
	sortedFileList = fileLengthSorter(smileFile);
	partDict = lineDivider(sortedFileList);
	global fileName;
	fileName = createFiles(partDict,smileFile,newPath);
    
    
if __name__ == "__main__":
	main(sys.argv[1]);

