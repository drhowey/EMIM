/************************************************************************
 * PREMIM, version 3.22
 * Copyright 2011-2016,
 * Richard Howey
 * Institute of Genetic Medicine, Newcastle University
 *
 * richard.howey@ncl.ac.uk
 * http://www.staff.ncl.ac.uk/richard.howey/
 *
 * This file is part of PREMIM, the pedigree file processing program for EMIM.
 *
 * PREMIM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PREMIM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PREMIM.  If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

/*! \file main.cpp
    \brief This file reads in the initial input files and options.
    
    This file also outputs usage instructions and program details.
*/

#include <iostream>
#include <ostream>
#include <sstream>
#include <cstdlib>
#include <time.h>

using namespace std; // initiates the "std" or "standard" namespace
 
#include "main.h"
#include "ProcessData.h"

bool outputToScreen = true; 

ofstream logFile;

//! Output program title to screen
void header()
{
	 out("\nPREMIM: Pedigree file processing program for EMIM, v3.22\n");
	 out("--------------------------------------------------------\n");
	 out("Copyright 2011-2016 Richard Howey, GNU General Public License, v3\n");
	 out("Institute of Genetic Medicine, Newcastle University\n\n");
};

//! Output program usage to screen
void usage()
{
	header();
		
#ifdef USING_GZIP
		out("Usage:\n\t ./premim [options] pedigreeFile [SNPNameFile]\n\n");
		out("  pedigreeFile -- a .ped, .bed or .gzip file\n");
#endif
#ifndef USING_GZIP
		out("Usage:\n\t premim [options] pedigreeFile [SNPNameFile]\n\n");
		out("  pedigreeFile -- a .ped or .bed file.\n");
#endif
		out("  SNPNameFile  -- a SNP description file (.map)\n\n");
		
		out("Data Processing Options:\n");
		out("  -a           -- estimate SNP allele frequencies (emimmarkers.dat)\n");
		out("  -s n dir     -- split SNP output into groups of size n in directory dir\n");
		out("  -n name      -- add a name to the end of the output files\n");
		out("  -xa          -- allow extra affected case trios per pedigree\n");
		out("  -xu          -- allow extra unaffected control matings per pedigree\n");
		out("  -pb file     -- add a proband file, listing affected subjects of interest\n");
		out("  -rmaj        -- use major allele as risk allele, \"2\"\n");
		out("  -rout file   -- output the risk alleles to file\n");
		out("  -rfile file  -- use alleles in file as risk alleles\n");
		out("  -log file    -- name of log file\n\n");

		out("Improved Imprinting using Haplotype Estimates:\n");
		out("  -ihap                   -- estimate haplotypes for imprinting\n");
		out("  -ihap-noadj             -- do not adjust estimated duo cell counts\n");
		out("  -ihap-miss-thres h      -- maximum missing data threshold, h, for trios and duos (default h=0.5)\n");
		out("  -shapeit com            -- shapeit command (default \"shapeit2\")\n");
		//out("  -shapeit-sams s         -- shapeit samples from posterior, 1 = best estimate\n");
		out("  -shapeit-thread t       -- shapeit thread option: --thread t (default t=12)\n");
		out("  -shapeit-mcmc-ops b p m -- shapeit options: --burn b --prune p --main m (default b=7, p=8, m=20)\n");
		out("  -shapeit-model-ops s w  -- shapeit options: --states s --window w (defualt s=100, w=2)\n");
		out("  -shapeit-ops \"options\"  -- other shapeit options, use \"\" to surround options\n\n");

		out("Parameter File Options:\n");
		out("  -cg          -- child genotype\n");
		out("  -ct          -- child trend (default)\n");
		out("  -mg          -- mother genotype\n");
		out("  -mt          -- mother trend\n");
		out("  -im          -- imprinting, maternal\n");
		out("  -ip          -- imprinting, paternal\n");
		out("  -imw         -- imprinting, maternal Wienberg\n");
		out("  -ipw         -- imprinting, paternal Wienberg\n\n");
		out("output is emimparams.dat and 11(or 12) other .dat files for input into EMIM.\n\n");
		out("Other Options:\n");
		out("  -so          -- suppress output to screen\n\n");
		out("Extra useful file manipulation operations (for parallel SNP analysis):\n");
		out("  premim -fm n markersFile.dat dir\n");
		out("  -- create marker files with n SNPs in dir using markersFile.dat\n\n");
		out("  premim -fr dir\n");
		out("  -- create results files emimresults.out (and emimsummary.out) from \n");
		out("     emimresults*.out (and emimsummary*.out) in dir\n\n");
};


//! Converts an integer to a string
string toString(int & i)
{
	ostringstream aStringStream;
	aStringStream << i;

	return aStringStream.str();
};

//! Returns a string of the run time
string getTime(const double & t)
{
	double time = t;
	int days = 0;
	int hours = 0;
	int minutes = 0;
	int seconds = 0;

	string ans = "";
	days = (int) (time / 86400); time -= days*86400;
	hours = (int) (time / 3600); time -= hours*3600;
	minutes = (int) (time / 60); time -= minutes*60;
	seconds = (int) time;

	if(days == 1) ans += "1 day";
	else if(days > 0) ans += toString(days) + " days";

	if(hours > 0)
	{
		if(days != 0)
		{
			if(minutes == 0 && seconds == 0) ans += " and ";
			else ans += ", ";
		};

		if(hours == 1) ans += "1 hour";
		else ans += toString(hours) + " hours";
	};

	if(minutes > 0)
	{
		if(ans != "")
		{
			if(seconds == 0) ans += " and ";
			else ans += ", ";
		};

		if(minutes == 1) ans += "1 minute";
		else ans += toString(minutes) + " minutes";
	};

	if(seconds > 0)
	{
		if(ans != "")
		{
			ans += " and ";			
		};

		if(seconds == 1) ans += "1 second";
		else ans += toString(seconds) + " seconds";
	};

	if(ans == "") ans = "less than one second";

	return ans;
};

//! The start of the program
int main(int argc, char * argv[])
{
	time_t start,end;
	double dif;
	time(&start);

	int argcount = 1;
	string fileName = "";
	string mapFileName = "";
	string logFileName = "";
	string probandFileName = "";
	string riskAlleleInFileName = "";
	string riskAlleleOutFileName = "";
	bool estimateAlleleFreq = false;
	bool extraAffectedTrios = false;
	bool extraUnaffectedTrios = false;
	bool childGenotype = false;
	bool childTrend = false;
	bool motherGenotype = false;
	bool motherTrend = false;
	bool imprintingMaternal = false;
	bool imprintingPaternal = false;
	bool imprintingMaternalWeinberg = false;
	bool imprintingPaternalWeinberg = false;
	bool splitMarkersFile = false;
	bool makeResultsFiles = false;
	bool useMajorAlleleAsRisk = false;
	bool pooHaps = false;
	bool adjustedCounts = true;
	string hapMCMCOptions = " --burn 7 --prune 8 --main 20"; 
	string hapModelOptions = " --states 100 --window 2"; 
	string threadOptions = " --thread 12";
	string otherShapeitOtherOptions = "";
	string shapeitCommand = "shapeit2";
	string temp;
	double missThres = 0.5;
	unsigned int noSamples = 1;

	outputToScreen = true;

	string option;
	string outputDirectory = "";
	string markersFile;
	string endName = "";
	unsigned int splitSNPOutput = 0;

	//set given options
	while(argcount < argc && argv[argcount][0] == '-')
    {
		option = argv[argcount];
		if(option == "-a") estimateAlleleFreq = true;	
		else if(option == "-xa") extraAffectedTrios = true;	
		else if(option == "-xu") extraUnaffectedTrios = true;	
		else if(option == "-cg") childGenotype = true;			
		else if(option == "-ct") childTrend = true;			
		else if(option == "-mg") motherGenotype = true;			
		else if(option == "-mt") motherTrend = true;			
		else if(option == "-im") imprintingMaternal = true;			
		else if(option == "-ip") imprintingPaternal = true;	
		else if(option == "-ihap") pooHaps = true;
		else if(option == "-ihap-noadj") {pooHaps = true; adjustedCounts = false;}
		else if(option == "-imw") imprintingMaternalWeinberg = true;			
		else if(option == "-ipw") imprintingPaternalWeinberg = true;
		else if(option == "-so") outputToScreen = false;
		else if(option == "-rmaj") useMajorAlleleAsRisk = true;
		else if(option == "-n")
		{
			argcount++; if(argcount >= argc) break;
			endName = argv[argcount];
		}
		else if(option ==  "-fm")
		{			
			argcount++; if(argcount >= argc) break;
			splitSNPOutput = atoi(argv[argcount]);
			argcount++; if(argcount >= argc) break;
			markersFile = argv[argcount];
			argcount++; if(argcount >= argc) break;
			outputDirectory = argv[argcount];
			splitMarkersFile = true;
		}
		else if(option ==  "-fr")
		{			
			argcount++; if(argcount >= argc) break;
			outputDirectory = argv[argcount];
			makeResultsFiles = true;
		}
		else if(option == "-s")
		{
			//get the number of SNPs to put into each file
			argcount++; if(argcount >= argc) break;
			splitSNPOutput = atoi(argv[argcount]);
			//get the name of directory to ouput the files into
			argcount++; if(argcount >= argc) break;
			outputDirectory = argv[argcount];
		}
		else if(option ==  "-pb")
		{			
			argcount++; if(argcount >= argc) break;
			probandFileName = argv[argcount];			
		}
		else if(option ==  "-rout")
		{			
			argcount++; if(argcount >= argc) break;
			riskAlleleOutFileName = argv[argcount];			
		}
		else if(option ==  "-rfile")
		{			
			argcount++; if(argcount >= argc) break;
			riskAlleleInFileName = argv[argcount];			
		}
		else if(option ==  "-log")
		{			
			argcount++; if(argcount >= argc) break;
			logFileName = argv[argcount];			
		}
		else if(option ==  "-ihap-miss-thres")
		{			
			argcount++; if(argcount >= argc) break;
			missThres = atof(argv[argcount]); //missing threshold			
		}
		else if(option ==  "-shapeit")
		{			
			argcount++; if(argcount >= argc) break;
			shapeitCommand = argv[argcount]; //shapeit command			
		}
		else if(option ==  "-shapeit-mcmc-ops")
		{						
			argcount++; if(argcount >= argc) break;
			temp = argv[argcount];
			hapMCMCOptions = " --burn " + temp; //set no of burn ins
			argcount++; if(argcount >= argc) break;
			temp = argv[argcount];
			hapMCMCOptions = hapMCMCOptions + " --prune " + temp; //set no of prunes
			argcount++; if(argcount >= argc) break;
			temp = argv[argcount];
			hapMCMCOptions = hapMCMCOptions + " --main " + temp; //set no of main iterations				
		}
		else if(option ==  "-shapeit-thread")
		{			
			argcount++; if(argcount >= argc) break;
			temp = argv[argcount];			
			threadOptions = " --thread "+ temp; //set no of thread				
		}
		else if(option ==  "-shapeit-model-ops")
		{		
			argcount++; if(argcount >= argc) break;
			temp = argv[argcount];			
			hapModelOptions = " --states "+ temp; //set states
			argcount++; if(argcount >= argc) break;
			temp = argv[argcount];
			hapModelOptions = hapModelOptions + " --window " + temp; //set window in kB							
		}
		else if(option == "-shapeit-sams")
		{
			argcount++; if(argcount >= argc) break;
			noSamples =  atoi(argv[argcount]);
		}
		else if(option ==  "-shapeit-ops")
		{	
			argcount++; if(argcount >= argc) break;
			otherShapeitOtherOptions = argv[argcount];			
		}
		else
		{
    		outErr("\nUnrecognised command line switch: "); outErr(option); outErr("\n");
			usage();
    		exit(1);
		};

		argcount++;
	};

	if(logFileName == "") logFileName = "premim" + endName + ".log";
	logFile.open(logFileName.c_str());

	//check that the given options are valid
	//allow two types of imprinting analysis if no child analysis is set
	if(childGenotype && childTrend)
	{
		outErr("\nOnly one type of child analysis is permitted!\n");
		exit(1);
	}
	else if(motherGenotype && motherTrend)
	{
		outErr("\nOnly one type of mother analysis is permitted!\n");
		exit(1);
	}
	else if(useMajorAlleleAsRisk && riskAlleleInFileName != "")
	{
		outErr("\nThe risk alleles may only be defined one way!\n");
		exit(1);
	}	
	else if(
		!((!childGenotype && !childTrend) && ((imprintingMaternal && imprintingPaternal)
		          || (imprintingMaternalWeinberg && imprintingPaternalWeinberg)))
		&&
		((imprintingMaternal && imprintingPaternal) || (imprintingMaternal && imprintingMaternalWeinberg)
		|| (imprintingMaternal && imprintingPaternalWeinberg) || (imprintingPaternal && imprintingMaternalWeinberg)
		|| (imprintingPaternal && imprintingPaternalWeinberg) || (imprintingMaternalWeinberg && imprintingPaternalWeinberg))		
		)
	{
		outErr("\nOnly one type of imprinting analysis is permitted (when performing a child analysis)!\n");
		exit(1);
	}
	else if(!childGenotype && !childTrend && !motherGenotype &&
		!motherTrend && !imprintingMaternal && !imprintingPaternal &&
		!imprintingMaternalWeinberg && !imprintingPaternalWeinberg && argc >= 2)
	{
		//set the default option if no options are set
		childTrend = true;
	};

	if(missThres < 0 || missThres >= 1)
	{
		outErr("\nThe maximum missing data threshold for parent-of-origin analysis must be between 0 and 1!\n");
		exit(1);
	}
	else if(extraAffectedTrios && pooHaps)
	{
		outErr("\nIt is not possible to use the extra trios option when using SHAPEIT2 to estimate the parent-of-origin of alleles.\n");
		outErr("The validity of parent-of-origin analyses using extra trios has not yet been investigated.\n");
		outErr("If you really want to do this (at your own risk), then manually create the extra trios in a .bed pedigree file and only use the ihap option.\n\n");
		exit(1);
	};
	
	if(argcount < argc || splitMarkersFile || makeResultsFiles)
	{
		header();

		out("Log file: "); out(logFileName); out("\n");

		if(!splitMarkersFile && !makeResultsFiles)
		{
			//get the file pedigree name
			fileName = argv[argcount++];

			//set the .map file if it is given
			if(argcount < argc) mapFileName = argv[argcount++];
		};

		//create an object to do the processing of the pedigree data
		ProcessData pData(extraAffectedTrios, extraUnaffectedTrios, childGenotype, childTrend, motherGenotype,
			motherTrend, imprintingMaternal, imprintingPaternal,
			imprintingMaternalWeinberg, imprintingPaternalWeinberg, estimateAlleleFreq, useMajorAlleleAsRisk, splitSNPOutput,
			outputDirectory, endName, probandFileName, riskAlleleInFileName, riskAlleleOutFileName, pooHaps, adjustedCounts, shapeitCommand, threadOptions, hapMCMCOptions, hapModelOptions, otherShapeitOtherOptions, noSamples, missThres);

		if(splitMarkersFile) pData.createMarkerFiles(splitSNPOutput, markersFile, outputDirectory, endName);
		else if(makeResultsFiles) pData.createResultsFiles(outputDirectory, endName);
		else
		{
			//process the pedigree data
			pData.process(fileName, mapFileName);
		};
	}
	else
	{
		usage();
	};

	time(&end);
	dif = difftime(end, start);
	out("Run time: "); out(getTime(dif)); out("\n\n");

	logFile.close();
};

