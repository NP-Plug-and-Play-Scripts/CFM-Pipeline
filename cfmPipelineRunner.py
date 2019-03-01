#!/usr/bin/env python3

"""
This script is the one you want to run if you have data from the NP-DB and want ran through CFM_predict so you get files ready for MS2LDA.
It calls all the other scripts and runs them after eachother.

IMPORTANT: its suggested to use the .sh instaler for CFM-ID found on this Github account. https://github.com/NP-Plug-and-Play-Scripts/Bash-scripts.git
This way the expected file structure is already in place. if not try to create a folder containing a folder with the following file structure.
(a tab means going in to a file)

-CFM_workplace (main folder for the workplace)
    -rdKit (library for cfm)
    -boost (library for cfm)
    -lpsolve (library for cfm)
    -cfm (is the cfm instalation and should contain a bin folder with runnable cfm tools).
        -bin (contains the runable cfm tools)
    -cfmData (is a folder that should contain 3 folders, and a cfm config file).
        -smileFile (location for the smiles data that you wish to run through this pipeline.)
        -results (location for the CFM_data results)
        -params_metab_ce_cfm  (contains the trained models for cfm_predict)
        -param_config.txt (contains the settings for cfm_predict)
    
also remember to put all the scripts below in the same folder or expect some issues.

Made by: Rutger Ozinga 
Last edit: 30/1/2019
"""

import os;
import re;
import sys;
import fileSplitter;
import inchiKeyCreator3;
import CFMrunner;
import spectraNormalizer;
import tandemMS_Merger;
import spectraDataEditor;
import fileMerger;
import smiles_neutralizer;
import smilesSplitter;
import optionFileParser;

def main(workplacePath,fileName,optionFile):
	exportString = "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:{0}boost_1_67_0/lib:{0}rdkit-Release_2016_03_1/lib:{0}lp_solve_5.5/lpsolve55/bin/ux64".format(workplacePath);
	os.environ["LD_LIBRARY_PATH"] = exportString;
	os.system(exportString);
	cfmPath = workplacePath + "cfm/bin/";
	cfmDataPath = workplacePath + "cfmData/";
	smilePath = cfmDataPath + "smileFile/";
	resultPath = cfmDataPath + "results/";
	"""Preparation"""
	#Obtain the options from the option file
	optionDict = optionFileParser.optionParser(optionFile);
	molConvertPath = optionDict["molconvert_path"];
	"""STEP 1"""
	#neutralize the smiles
	print("neutralizing Smiles");
	neutralizedPath = smiles_neutralizer.main(smilePath + fileName);
	neutralCSV = neutralizedPath.split("/")[-1];
	print("Smiles Neutralized!")
	#split files
	"""STEP 2"""
	fileSplitter.main(smilePath + neutralCSV);
	print("File splitting done!");
	fileName = fileSplitter.fileName;
	"""STEP 3"""
	#create inchi keys and add them to a file containing a DB id and smiles string
	print("creating inchiKeys");
	inchiKeyCreator3.main(molConvertPath, smilePath, neutralCSV);
	print("inchiKeys created!");
	"""STEP 4"""
	#split the Smiles if they have multiple molecules in the smile.
	print("starting with the splitting");
	smilesSplitter.main(smilePath, fileName);
	print("smiles split and longest ones saved");
	"""STEP 5"""
	#run CFM_ID
	print("getting ready to run CFM-ID");
	CFMrunner.main(cfmPath,cfmDataPath,fileName,resultPath,optionDict);
	print("Done running CFM_ID!");
	"""STEP 6"""
	#run spectraNormalizer
	print("normalizing spectra");
	spectraNormalizer.main(resultPath,fileName);
	print("done normalizing");
	"""STEP 7"""
	#run spectraMerger
	print("merging the spectra")
	tandemMS_Merger.main(resultPath,fileName, optionDict);
	print("merging done!");
	"""STEP 8"""
	#run spectraDataEditor
	print("adding extra data to the spectra");
	spectraDataEditor.main(smilePath,resultPath,fileName);
	print("Complete the data is now ready for MS2LDA");
	"""STEP 9"""
	#put all files together again.
	print("combining files");
	fileMerger.main(resultPath,fileName);
	print("Files now combined!");
	print("End of pipeline, have a nice day!");

if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2], sys.argv[3]);
