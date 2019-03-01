#!/usr/bin/env python
import sys;
"""
parses a option file in to a dict. Starts parsing after Options_Start has passed and stops when Options_End comes by
"""
def optionParser(optionPath):
	optionDict = {};
	parseOptions = False;
	for optionLine in open(optionPath):
		optionLine = optionLine.strip();
		if optionLine == "Options_Begin":
			parseOptions = True;
		elif optionLine == "Options_End":
			parseOptions = False;
		#if parseOptions is true start splitting and adding lines to the optionDict.
		if "=" in optionLine:
			if parseOptions:
				print(optionLine);
				splittedOption = optionLine.split("=");
				optionDict[splittedOption[0]] = splittedOption[1];
	return optionDict;

if __name__ == "__main__":
	optionParser(sys.argv[1]);
