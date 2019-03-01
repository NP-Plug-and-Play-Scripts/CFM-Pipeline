#!/usr/bin/env python3
"""
Runs all the files in a directory with a given name through cfm_predict. required a path to the cfm-id bin, results folder for output
cfmData folder for input and the name of the file (name is taken from the file in the fileSplitter.py script normaly looks like ether_131_part_)
then runs all of the files that correspond to the pattern on seperate cores in cfm_predict (10 by default will add the option to change it later).

CFM_pipeline Step 4: 
Third step in the pipeline (although it can be second too since its only dependant on fileSplitter.py). Takes the NP-ID smile combinations 
and puts them in cfm_predict to genereate the in-silico spectra. In time i hope to add more options to change the settings more easily. next
next script is spectraNormalizer.

Made by: Rutger Ozinga 
Last edit: 10/10/2018
"""
import os;
import sys;
import multiprocessing;
import re;

"""
Runs CFM ID with the given info
cfmPath = path to the cfm-id runable folder (../cfm/bin)
cfmData = path to the cfmData folder
inName = the name of the Input file
outPath = path to which the output is written
cfmOptionDict = a dict containing the options for the pipeline.
"""
def runCFM(cfmPath,cfmData,inName,outPath,optionDict):
	#the config file this contains the options for cfm.
	config = optionDict["config"];
	#the parameters for the likelyhood of breaking this contains the pretrained model.
	parameters = optionDict["parameters"];
	#probability threshold 
	probThresh = optionDict["prob_thresh"];
	cfmMode = cfmPath + "cfm-predict";
	cfmParam = cfmData + parameters;
	cfmConfig = cfmData + config;
	filePath = cfmData + "smileFile/" + inName;
	command = "{0} {1} {2} {3} {4} 0 {5}".format(cfmMode,filePath,probThresh,cfmParam,cfmConfig,outPath);
	os.system(command);
	#cfmPath/cfm-predict ${cfmData}${smileDir}${newFileName}${value}.txt 0.001 $cfmData/params_metab_ce_cfm/param_output0.log $cfmData/param_config.txt 0 $cfmData/results/${newFileName}${value}_output.mgf&

"""
get all the files in a given file path that corresponds to a given pattern.
cfmData = path to the cfmData file
filePattern = patern that matches to the first part of the file
"""
def getFileList(cfmData,filePattern):
	pattern = re.escape(filePattern) + r'[0-9]{2}.txt';
	fileList = [f for f in os.listdir(cfmData+"/smileFile/") if re.search(pattern,f)];
	#sort it so the files go from 00 to 09;
	fileList.sort();
	return fileList;


"""
Main method runs all the functions.
Runs all the differe file parts in seperate jobs and waits for them to all be finished before moving on.
cfmPath = path to the cfm-d bin folder
cfmData = path to the cfmData folder
fileName = pattern that matches the first part of the file "ethers_131_part_" or something similair.
outputPath = path to the results folder.
"""
def main(cfmPath, cfmData,fileName,outputPath,optionDict):
	fileList = getFileList(cfmData,fileName);
	#list of running jobs
	jobs = []
	for x in range(len(fileList)):
		newFileName = fileList[x].replace(".txt","_output.mgf");
		newFilePath = outputPath + newFileName;
		#creates a process object out of runCFM ands adds the job to a list of jobs. this way 
		job = multiprocessing.Process(target=runCFM, args=(cfmPath,cfmData,fileList[x],newFilePath,optionDict))
		jobs.append(job);
		#start the job.
		job.start();
	#This loop makes it so the script doesnt continue before all processes are finished.
	for job in jobs:
		job.join();

if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]);

