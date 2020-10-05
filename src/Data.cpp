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


/*! \file Data.cpp
    \brief This file contains functions for processing basic genotype infomation.
       
*/

#include "Data.h"
#include "main.h"
#include "ModelFitting.h"

#include <string>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <cstdlib>

using namespace std; // initiates the "std" or "standard" namespace

//! Get the ID of a subject, family ID + subject ID maps to an ordinal number over all subjects. 
unsigned int MapIds::getId(const string & name)
{
  map<string, unsigned int>::const_iterator i = idMap.find(name);
  if(i != idMap.end()) return i->second;

  //if id does not exist create an id
  unsigned int id = (unsigned)idMap.size() + 1;
  idMap[name] = id;
  
  return id;
};

//! Get the name of a subject, given the id.
string MapIds::getName(const unsigned int & id)
{
	//loop thro' names to find the id, this is not too quick but is only used for reporting errors/warnings
	for(map<string, unsigned int>::const_iterator i = idMap.begin(); i != idMap.end(); ++i)
	{
		if(id == i->second) return i->first;
	};

	return "ID not found";
};

//! Get the allele ID (0 or 1) from the allele name, e.g. for a SNP: (A, G) <-> (0, 1).
bool MapSnpAlleleIds::getSnpAlleleId(const unsigned int & snpId, string & alleleName)
{

	map<unsigned int, map<string, bool> >::iterator i = mapSnpAlleles.find(snpId); 

	//if allele name to id (0 or 1) exists then return it else create id
	if(i != mapSnpAlleles.end())
	{
		map<string, bool>::iterator a = i->second.find(alleleName);
		if(a != i->second.end())
		{
			return a->second;
		}
		else
		{
			//the other allele value is already set to 0 so set this allele to 1
			i->second[alleleName] = 1;

			//check for non-SNPs
			if(i->second.size() > 2)
			{				
				outErr("Warning: SNP "); outErr(snpId); outErr(" has "); outErr(i->second.size()); outErr(" allele values and is not a SNP!\n");				
			};

			return 1;
		};

	};
	
	//set the first allele name encountered to 0
	map<string, bool> alleleIds;
	alleleIds[alleleName] = 0;

	mapSnpAlleles[snpId] = alleleIds;
	
	return 0;
};

//! Get the allele names for a SNP.
map<string, bool> MapSnpAlleleIds::getSnpAlleleNames(const unsigned int & snpId) const
{
	map<unsigned int, map<string, bool> >::const_iterator msa = mapSnpAlleles.find(snpId);

	if(msa != mapSnpAlleles.end()) return msa->second;
	else
	{
		map<string, bool> dummy;
		dummy["NoData "] = 0;
		dummy["NoData"] = 1;
		return dummy;
	};
};

//! Adds genotype data and always stores 11, 12, 21 and 22 as 00, 01 and 11
void GenotypeData::addGenotype(const unsigned int & snpId, const bool & allele1, const bool & allele2)
 {
   bool a1 = allele1;
   bool a2 = allele2;
   if(allele1 == 1 && allele2 == 0)
     {
       a1 = 0; a2 = 1; //ensure that the genotype data is always stored with the first allele first
     };
   data[snpId] = make_pair(a1, a2);
 };

//! Finds genotype infomation for a given SNP.
pair<bool, bool> GenotypeData::getGenotype(const unsigned int & snpId) const
{
	map<unsigned int, pair<bool, bool> >::const_iterator g = data.find(snpId);

	if(g != data.end())
	{
		return g->second;
	}
	else
	{		
		//unknown genotype
		//error
		outErr("Genotype not found for SNP "); outErr(snpId); outErr("!\n");
		exit(1);
	};
};

//! Checks if genotype data for a given SNP exists or not.
bool GenotypeData::genotypeDataPresent(const unsigned int & snpId) const
{
  map<unsigned int, pair<bool, bool> >::const_iterator g = data.find(snpId);

  if(g != data.end()) return true;

  return false;
};

//! Checks if the SNP that has been labelled allele "1" is the minor allele or not.
bool CountAlleles::isSnpAllele1Minor(const unsigned int & snpId) const
{
	map<unsigned int, bool>::const_iterator i = isAllele1Minor.find(snpId);

	if(i != isAllele1Minor.end()) return i->second;
	else return false; //case where there is no genotype infomation
};

//! Checks if the alleles were switched to match risk allele requirements.
bool CountAlleles::getWereAllelesReversed(const unsigned int & snpId) const
{
	map<unsigned int, bool>::const_iterator i = wereAllelesReversed.find(snpId);

	if(i != wereAllelesReversed.end()) return i->second;
	else return false; //case where there is no genotype infomation
};

//! Checks if there is any SNP data.
bool CountAlleles::isSNPDataMissing(unsigned int & snpId) const
{
	map<unsigned int, AlleleCounts>::const_iterator i = allele1Counts.find(snpId);

	if(i != allele1Counts.end()) return ((i->second.totalSubjectsAffected + i->second.totalSubjectsUnaffected) == 0);
	return true;
};

//! Set whether allele "1" is minor for each SNP or not
void CountAlleles::updateIsAllele1Minor()
{
	//check if the allele count is less than half of the maximum possible allele count (= (2 * no of Subjects)/2 = no of Subjects)
	for(map<unsigned int, AlleleCounts>::iterator ac = allele1Counts.begin(); ac != allele1Counts.end(); ++ac)
	{		
		isAllele1Minor[ac->first] = (ac->second.unaffectedAllele1 + ac->second.affectedAllele1) < (ac->second.totalSubjectsUnaffected + ac->second.totalSubjectsAffected);
	};
};

//! Set whether allele "1" is minor for each SNP or not
void CountAlleles::updateIsAllele1Minor(unsigned int & snpID)
{
	map<unsigned int, AlleleCounts>::const_iterator ac = allele1Counts.find(snpID);

	//check if the allele count is less than half of the maximum possible allele count (= (2 * no of Subjects)/2 = no of Subjects)
	if(ac != allele1Counts.end())
	{		
		isAllele1Minor[snpID] = (ac->second.unaffectedAllele1 + ac->second.affectedAllele1) < (ac->second.totalSubjectsUnaffected + ac->second.totalSubjectsAffected);
	};
};

//! Adds to the allele "1" count taking into account affected/unaffected subjects and number of defined genotypes
void CountAlleles::addToAllele1Count(const unsigned int & snpId, const unsigned int & i, const bool & affected)
{
	map<unsigned int, AlleleCounts>::iterator ac = allele1Counts.find(snpId);
	
	//if the allele count object does not exist then it is create else the existing total is added to
	//we need to account for the total since some genotypes may be missing
	if(ac == allele1Counts.end())
	{			
		if(affected)
		{
			AlleleCounts alleleCounts(i, 0, 1, 0); 
			allele1Counts[snpId] = alleleCounts;
		}
		else
		{
			AlleleCounts alleleCounts(0, i, 0, 1); 
			allele1Counts[snpId] = alleleCounts;
		};	
	}
	else
	{
		if(affected)
		{
			ac->second.affectedAllele1 += i;
			ac->second.totalSubjectsAffected++;			
		}
		else
		{
			ac->second.unaffectedAllele1 += i;
			ac->second.totalSubjectsUnaffected++;
		};	
	};
	
};

//! If the SNP info is to be split across different files then check to see if a new file is needed.
void updateMarkersFile(ofstream & markersFile, const unsigned int & snpId, const unsigned int & splitSNPOutput, const string & outputDirectory, const string & endName)
{	
	if(splitSNPOutput > 0)
	{
		if((snpId-1)%splitSNPOutput == 0)
		{
			if(markersFile.is_open()) markersFile.close();			
			unsigned int fileNo = (unsigned int)((snpId-1)/splitSNPOutput) + 1;
			ostringstream nextFilename;
			nextFilename << outputDirectory << "emimmarkers"<< endName << fileNo << ".dat";
			string fname = nextFilename.str();
			markersFile.open(fname.c_str());
		};
	};
};

//! Writes the allele frequency estimates to file.
void CountAlleles::outputMarkersFile(ofstream & markersFile, const unsigned int & noOfSnps, const unsigned int & splitSNPOutput, const string & outputDirectory, const string & endName,
	bool & useMajorAlleleAsRisk, map<unsigned int, string> & riskAlleleNames, MapSnpAlleleIds & snpAlleleIds, const bool & useGivenRiskAlleles)
{	
	double freqEstimate = 0.01;
	bool allele1IsRiskAllele;
	map<unsigned int, string>::const_iterator ra;
	map<string, bool> snpAls;
	map<string, bool>::const_iterator al;

	map<unsigned int, AlleleCounts>::const_iterator ac = allele1Counts.begin();
	for(unsigned int snpId = 1; snpId <= noOfSnps; ++snpId)
	{
		freqEstimate = 0.01;
		if(ac != allele1Counts.end() && snpId == ac->first)
		{
			if((ac->second.totalSubjectsUnaffected + ac->second.totalSubjectsAffected) != 0) freqEstimate = 0.5*(((double)(ac->second.unaffectedAllele1 + ac->second.affectedAllele1))/((double)(ac->second.totalSubjectsUnaffected + ac->second.totalSubjectsAffected))); 

			++ac;
		};

		//set using given risk alleles, set by user file, -rfile
		if(useGivenRiskAlleles)
		{
			ra = riskAlleleNames.find(snpId); //snpNo, allele name

			if(ra == riskAlleleNames.end())
			{
				out("Risk allele not found for SNP number "); out(snpId); out("!\n");
				exit(1);
			};

			snpAls = snpAlleleIds.getSnpAlleleNames(snpId); //map<unsigned int, map<string, bool> > //SNP ID no., <allele name,  allele id - 0 or 1>
			al = snpAls.begin(); //consider the first allele

			//check if risk allele is same as the first allele and this allele is allele 1 (for PREMIM)
			if(ra->second == al->first && !al->second) allele1IsRiskAllele = true;
			else
			{
				//check if second allele is allele 1 (in PREMIM) and is risk allele
				al++;
				if(al != snpAls.end() && ra->second == al->first && !al->second) allele1IsRiskAllele = true;
				else allele1IsRiskAllele = false;
			};
			
		}
		else
		{

			//Check if the first allele is the risk allele. First allele according to PREMIM, i.e. the first allele encounted in file this may or may not be allele "1" as given in EMIM docs
			allele1IsRiskAllele = (isSnpAllele1Minor(snpId) && !useMajorAlleleAsRisk) || (!isSnpAllele1Minor(snpId) && useMajorAlleleAsRisk);

		};

		if(!allele1IsRiskAllele) freqEstimate = 1.0 - freqEstimate;

		if(freqEstimate < 0.01) freqEstimate = 0.01;
		if(freqEstimate > 0.99) freqEstimate = 0.99;

		updateMarkersFile(markersFile, snpId, splitSNPOutput, outputDirectory, endName);
		markersFile << snpId << "\t" << freqEstimate <<"\n";
	};
	
	markersFile.close();
};

//! Get minor allele freqs using all data.
map<unsigned int, double> CountAlleles::getMinorAlleleFreqs(unsigned int & noOfSNPs)
{	
	map<unsigned int, double> freqs;
	double freq; //of allele1
	map<unsigned int, AlleleCounts>::const_iterator ac = allele1Counts.begin();
	for(unsigned int snpId = 1; snpId <= noOfSNPs; ++snpId)
	{
		
		if(ac != allele1Counts.end() && snpId == ac->first)
		{
			freq = 0.5*((double)(ac->second.unaffectedAllele1 + ac->second.affectedAllele1))
				             /((double)(ac->second.totalSubjectsUnaffected + ac->second.totalSubjectsAffected));		
			++ac;
		}
		else freq = 0;

		if(freq > 0.5) freq = 1 - freq;	
		freqs[snpId] = freq;
	};

	return freqs;
};

//! Outputs the estimated allele frequency for one SNP.
void CountAlleles::outputMarkersFileOneLine(ofstream & markersFile, const unsigned int & snpId, const unsigned int & splitSNPOutput, const string & outputDirectory, const string & endName,
	bool & useMajorAlleleAsRisk, map<unsigned int, string> & riskAlleleNames, MapSnpAlleleIds & snpAlleleIds, const bool & useGivenRiskAlleles)
{	
	//default frequency for very low estimates or SNPs with no estimate
	double freqEstimate = 0.01;
	bool allele1IsRiskAllele;
	map<unsigned int, string>::const_iterator ra;
	map<string, bool> snpAls;
	map<string, bool>::const_iterator al;

	//estimate the frequency using the unaffected subjects if possible otherwise use the affected subjects
	map<unsigned int, AlleleCounts>::const_iterator ac = allele1Counts.begin();

	if(ac != allele1Counts.end())
	{
		if((ac->second.totalSubjectsUnaffected + ac->second.totalSubjectsAffected) != 0) freqEstimate = 0.5*(((double)(ac->second.unaffectedAllele1 + ac->second.affectedAllele1))/((double)(ac->second.totalSubjectsUnaffected + ac->second.totalSubjectsAffected))); 
		else freqEstimate = 0;
	};
	
	//set using given risk alleles, set by user file, -rfile
	if(useGivenRiskAlleles)
	{
		ra = riskAlleleNames.find(snpId);

		if(ra == riskAlleleNames.end())
		{
			out("Risk allele not found for SNP number "); out(snpId); out("!\n");
			exit(1);
		};

		snpAls = snpAlleleIds.getSnpAlleleNames(snpId);			
		al = snpAls.begin(); //consider the first allele

		//check if risk allele is same as the first allele and this allele is allele 1 (for PREMIM)
		if(ra->second == al->first && !al->second) allele1IsRiskAllele = true;
		else
		{
			//check if second allele is allele 1 (in PREMIM) and is risk allele
			al++;
			if(al != snpAls.end() && ra->second == al->first && !al->second) allele1IsRiskAllele = true;
			else allele1IsRiskAllele = false;
		};	

	}
	else
	{
		//Check if the first allele is the risk allele. First allele according to PREMIM, i.e. the first allele encounted in file this may or may not be allele "1" as given in EMIM docs
		allele1IsRiskAllele = (isSnpAllele1Minor(snpId) && !useMajorAlleleAsRisk) || (!isSnpAllele1Minor(snpId) && useMajorAlleleAsRisk);
	};

	if(!allele1IsRiskAllele) freqEstimate = 1.0 - freqEstimate;

	if(freqEstimate < 0.01) freqEstimate = 0.01;
	if(freqEstimate > 0.99) freqEstimate = 0.99;

	updateMarkersFile(markersFile, snpId, splitSNPOutput, outputDirectory, endName);

	//write the frequency estimate to the file
	markersFile << snpId << "\t" << freqEstimate <<"\n";	
};

//! Adds to the count of different genotypes for case parent trios or case mother duos etc. etc.
void CountedGenotypeData::addItem(const unsigned int & snpId, unsigned int & groupId)
{

	map<unsigned int, map<unsigned int, unsigned int> >::iterator s = countedData.find(snpId);

	if(s != countedData.end())
	{
		map<unsigned int, unsigned int>::iterator g = s->second.find(groupId);
		if(g != s->second.end())
		{
			countedData[snpId][groupId]++;
		}
		else
		{
			countedData[snpId][groupId] = 1;
		};
	}
	else
	{
		map<unsigned int, unsigned int> groupCounts;
		groupCounts[groupId] = 1;
		countedData[snpId] = groupCounts;
	};

};


//! Adds to the count of different genotypes for case parent trios or case mother duos etc. etc.
void CountedGenotypeData::addItemPhased(const unsigned int & snpId, unsigned int & groupId)
{

	map<unsigned int, map<unsigned int, unsigned int> >::iterator s = countedDataPhased.find(snpId);

	if(s != countedDataPhased.end())
	{
		map<unsigned int, unsigned int>::iterator g = s->second.find(groupId);
		if(g != s->second.end())
		{
			countedDataPhased[snpId][groupId]++;
		}
		else
		{
			countedDataPhased[snpId][groupId] = 1;
		};
	}
	else
	{
		map<unsigned int, unsigned int> groupCounts;
		groupCounts[groupId] = 1;
		countedDataPhased[snpId] = groupCounts;
	};
};

//! Writes to file the counted genotype data
void CountedGenotypeData::writeSnpGroups(ofstream & fileOutPut, const unsigned int & snpId, map<unsigned int, map<unsigned int, unsigned int> >::const_iterator & snp, map<unsigned int, map<unsigned int, unsigned int> >::const_iterator & snpPh, const bool & oneSnp)
{
	if(totalSamplesPhased != 0)
	{
			writeSnpGroupsPhased(fileOutPut, snpId, snp, snpPh, oneSnp);
			return;
	};

	//output the SNP Id at the begining of each line
	fileOutPut << snpId << "\t";

	if(snp != countedData.end() && (oneSnp || snp->first == snpId))
	{
		map<unsigned int, unsigned int>::const_iterator grp = snp->second.begin();
		if(grp->first == 0)
		{
			undefinedSnpGroups = grp->second;
			++grp;
		};

		for(unsigned int grpNo = 1; grpNo <= noOfGroups; )
		{
			if(grp != snp->second.end() && grp->first == grpNo)
			{
				fileOutPut << grp->second;

				totalSnpGroups += grp->second;
				++grp;
			}
			else
				fileOutPut << "0";
			
			grpNo++;
			if(grpNo <= noOfGroups ) fileOutPut << " ";
		};

		fileOutPut << "\n";

		++snp;
	}
	else
	{
		for(unsigned int grpNo0 = 1; grpNo0 <= noOfGroups; )
		{
			fileOutPut << "0";
			grpNo0++;
			if(grpNo0 <= noOfGroups ) fileOutPut << " ";
		};
		fileOutPut << "\n";
	};

};

//! Writes to file the counted genotype data
void CountedGenotypeData::writeSnpGroupsPhased(ofstream & fileOutPut, const unsigned int & snpId, map<unsigned int, map<unsigned int, unsigned int> >::const_iterator & snp, map<unsigned int, map<unsigned int, unsigned int> >::const_iterator & snpPh, const bool & oneSnp)
{
	map<unsigned int, double> counts; //grpNo, count
	map<unsigned int, unsigned int>::const_iterator grp;
	
	unsigned int totalGroupsThisSNP = 0;
	unsigned int totalGroupsPhasedThisSNP = 0;

	//output the SNP Id at the begining of each line
	fileOutPut << snpId <<"\t";

	if(snp != countedData.end() && (oneSnp || snp->first == snpId))
	{
		grp = snp->second.begin();
		if(grp->first == 0)
		{
			undefinedSnpGroups = grp->second;
			++grp;
		};

		for(unsigned int grpNo = 1; grpNo <= noOfGroups; )
		{
			if(grp != snp->second.end() && grp->first == grpNo)
			{			
				counts[grpNo] += grp->second;
				totalSnpGroups += grp->second;
				totalGroupsThisSNP += grp->second;
				++grp;
			}
			else
				counts[grpNo] = 0;			
			
			grpNo++;
			if(grpNo <= noOfGroups ) fileOutPut << " ";
		};		

		++snp;
	}
	else
	{
		for(unsigned int grpNo0 = 1; grpNo0 <= noOfGroups; )
		{			
			counts[grpNo0] = 0;
			grpNo0++;		
		};
		
	};

	unsigned int totalSnpGroupsPhased = 0;
	double totalSamples = (double)(totalSamplesPhased);
	double phasedCount;
	map<unsigned int, double> adjustedPhasedCounts;
	map<unsigned int, double>::const_iterator ac;

	if(adjustedCounts)
	{
		map<unsigned int, map<unsigned int, double> >::const_iterator cdpa = countedDataPhasedAdjusted.find(snpId);
		if(cdpa != countedDataPhasedAdjusted.end())	adjustedPhasedCounts = cdpa->second;
	};

	//add phased counts
	if(snpPh != countedDataPhased.end() && (oneSnp || snpPh->first == snpId))
	{
		grp = snpPh->second.begin();
		if(grp->first == 0)
		{
			undefinedSnpGroups += grp->second;
			++grp;
		};

		for(unsigned int grpNo = 1; grpNo <= noOfGroups; )
		{
			if(grp != snpPh->second.end() && grp->first == grpNo)
			{	
				//set to adjusted count if it exists
				if(adjustedCounts)
				{
					ac = adjustedPhasedCounts.find(grpNo);
					if(ac != adjustedPhasedCounts.end()) phasedCount = ac->second/totalSamples;
					else phasedCount = (double)(grp->second)/totalSamples; 
				}
				else
					phasedCount = (double)(grp->second)/totalSamples;

				counts[grpNo] += phasedCount;
				totalSnpGroupsPhased += grp->second; //no need to adjust total as these will be the same
				totalGroupsPhasedThisSNP += grp->second;
				++grp;
			};			
			
			grpNo++;			
		};		

		++snpPh;
	};
	
	totalSnpGroups += totalSnpGroupsPhased/totalSamplesPhased; //should divide

	//output counts
	for(map<unsigned int, double>::const_iterator cts = counts.begin(); cts != counts.end(); )
	{
		fileOutPut << cts->second;
		++cts;
		fileOutPut << " ";
	};

	//output totals for normal EMIM and EMIM using shapeit2
	fileOutPut << totalGroupsThisSNP << " " << totalGroupsPhasedThisSNP << "\n";

};

//! Outputs results to file for one SNP.
void CountedGenotypeData::outputResultFileOneSnp(const unsigned int & snpId, const unsigned int & splitSNPOutput, const bool & pooHaps)
{
	//open a new file if splitting the output over a number of files and it is due
	if(splitSNPOutput > 0)
	{
		if(snpId == 1)
		{
			string fname0 = filename + "1.dat";
			fileOutput.open(fname0.c_str());
			if(!pooHaps) fileOutput << "snp\tcellcount 1-"<<noOfGroups<<"\n";
			else fileOutput << "snp\tcellcount 1-"<<noOfGroups<<" (+2 phased, +2 totals)\n";
		}
		else if((snpId-1)%splitSNPOutput == 0)
		{
			if(fileOutput.is_open()) fileOutput.close();
			unsigned int fileNo = (unsigned int)((snpId-1)/splitSNPOutput) + 1;
			ostringstream nextFilename;
			nextFilename << filename << fileNo << ".dat";
			string fname = nextFilename.str();
			fileOutput.open(fname.c_str());
			if(!pooHaps) fileOutput << "snp\tcellcount 1-"<<noOfGroups<<"\n";
			else fileOutput << "snp\tcellcount 1-"<<noOfGroups<<" (+2 phased, +2 totals)\n";
		};

	};

	map<unsigned int, map<unsigned int, unsigned int> >::const_iterator snp = countedData.begin();
	map<unsigned int, map<unsigned int, unsigned int> >::const_iterator snpPh = countedDataPhased.begin();

	writeSnpGroups(fileOutput, snpId, snp, snpPh, true); 

	//remember if something was found to write to parameter file
	if(hasData()) somethingCounted = true;
	countedData.clear(); //clear data for next SNP
};

//! Outputs a summary of counted genotype data.
void outputSummary(const string & groupName, string & filename, const unsigned int & noOfSnps, const unsigned int & totalNoOfPedigrees, const unsigned int & noOfGroups, const unsigned int & totalSnpGroups, const unsigned int & undefinedSnpGroups)
{	
	double avePerSnp = (((double)(totalSnpGroups))/(double)(noOfSnps));;

	out("File name: "); out(filename);			
	out("\nNumber of counted "); out(groupName); out(" (all SNPs): "); out(totalSnpGroups);
	out("\nAverage number of counted "); out(groupName); out(" (per SNP): "); out(avePerSnp); out("\n");
	if(noOfGroups > 3) {out("Number of uncounted (Mendelian error) "); out(groupName); out(": "); out(undefinedSnpGroups); out("\n\n");}
	else out("\n");
};

//! Outputs a result file for counted genotypes.
void CountedGenotypeData::outputResultsFile(const unsigned int & noOfSnps, const unsigned int & totalNoOfPedigrees, const unsigned int & splitSNPOutput)
{
	string snpFilename;
	
	//output all counted SNPs to one file
	if(splitSNPOutput == 0)
	{
		//snpFilename = filename+".dat";
		//fileOutput.open(snpFilename.c_str());
		//fileOutput << "snp\tcellcount 1-"<<noOfGroups<<"\n";
	
		map<unsigned int, map<unsigned int, unsigned int> >::const_iterator snp = countedData.begin();
		map<unsigned int, map<unsigned int, unsigned int> >::const_iterator snpPh = countedDataPhased.begin();
		
		for(unsigned int snpId = 1; snpId <= noOfSnps; ++snpId)
		{
			writeSnpGroups(fileOutput, snpId, snp, snpPh, false); 
		};

		fileOutput.close();
	}
	else
	{
		//output counted SNPs to different file
		unsigned int fileNo = 1;
		unsigned int snpId = 1;
		ostringstream snfn;

		map<unsigned int, map<unsigned int, unsigned int> >::const_iterator snp = countedData.begin();
		map<unsigned int, map<unsigned int, unsigned int> >::const_iterator snpPh = countedDataPhased.begin();

		do{
			snfn.str("");
			snfn << filename << fileNo << ".dat";
			snpFilename = snfn.str();

			fileOutput.open(snpFilename.c_str());
			fileOutput << "snp\tcellcount 1-"<<noOfGroups<<"\n";
		
			//add splitSNPOutput number of SNPs to the file, or less if there are not enough SNPs left
			for(unsigned int count = 1; (count <= splitSNPOutput && snpId <= noOfSnps); ++count)
			{
				writeSnpGroups(fileOutput, snpId, snp, snpPh, false);
				++snpId; 
			};			

			fileOutput.close();
			++fileNo;
		}while(snpId <= noOfSnps);
	};

	outputSummary(groupName, filename, noOfSnps, totalNoOfPedigrees, noOfGroups, totalSnpGroups, undefinedSnpGroups);
};

//! Outputs a summary of genotype data to file. 
void CountedGenotypeData::outputFileSummary(const unsigned int & noOfSnps, const unsigned int & totalNoOfPedigrees)
{	
	outputSummary(groupName, filename, noOfSnps, totalNoOfPedigrees, noOfGroups, totalSnpGroups, undefinedSnpGroups);
};

//! Outputs details of a subject to screen.
void Subject::outputDetails()
{
	out("Subject id: "); out(id); out("\n");
	out("Father id: "); out(fatherId); out("\n");
	out("Mother id: "); out(motherId); out("\n");
	out("Sex id: "); out(sex); out("\n");
	out("Affected id: "); out(affected); out("\n\n");
};

//! Returns the genotype group (1-15) a child parents trio belongs, except tweaked for phased data.
unsigned int CountedTrios::getGroupPhased(const unsigned int & snpId, bool & fatherAllele1, bool & fatherAllele2, bool & motherAllele1, bool & motherAllele2, bool & childAllele1, bool & childAllele2)
{
	if(motherAllele1 == 1 && motherAllele2 == 1)
	{
		if(fatherAllele1 == 1 && fatherAllele2 == 1)
		{
			if(childAllele1 == 1 && childAllele2 == 1)
			{
				return 1;
			};
			//else error
		}
		else if(fatherAllele1 == 0 && fatherAllele2 == 0)
		{
			if((childAllele1 == 0 && childAllele2 == 1) || (childAllele1 == 1 && childAllele2 == 0))
			{
				return 6;
			};
			//else error
		}
		else //father is 0 / 1
		{
			if((childAllele1 == 0 && childAllele2 == 1) || (childAllele1 == 1 && childAllele2 == 0))
			{
				return 3;
			}
			else if(childAllele1 == 1 && childAllele2 == 1)
			{
				return 2;
			};
			//else error
		};
    }
	else if(motherAllele1 == 0 && motherAllele2 == 0)
	{
		if(fatherAllele1 == 1 && fatherAllele2 == 1)
		{
			if((childAllele1 == 0 && childAllele2 == 1) || (childAllele1 == 1 && childAllele2 == 0))
			{
				return 7;
			};
			//else error
		}
		else if(fatherAllele1 == 0 && fatherAllele2 == 0)
		{

			if(childAllele1 == 0 && childAllele2 == 0)
			{
				return 15;
			};
			//else error
		}
		else//father is 0 / 1
		{
			if((childAllele1 == 0 && childAllele2 == 1) || (childAllele1 == 1 && childAllele2 == 0))
			{
				return 13;
			}
			else if(childAllele1 == 0 && childAllele2 == 0)
			{
				return 14;
			}
			//else error
		};
	}
	else //mother is 0 / 1
	{
		if(fatherAllele1 == 1 && fatherAllele2 == 1)
		{
			if((childAllele1 == 0 && childAllele2 == 1) || (childAllele1 == 1 && childAllele2 == 0))
			{
				return 5;
			}
			else if(childAllele1 == 1 && childAllele2 == 1)
			{
				return 4;
			};
			//else error
		}
		else if(fatherAllele1 == 0 && fatherAllele2 == 0)
		{
			if((childAllele1 == 0 && childAllele2 == 1) || (childAllele1 == 1 && childAllele2 == 0))
			{
				return 11;
			}
			else if(childAllele1 == 0 && childAllele2 == 0)
			{
				return 12;
			};
			//else error
		}
		else//father is 0 / 1
		{
			if(childAllele1 == 0 && childAllele2 == 1)
			{						
				return 17; 
			}
			else if(childAllele1 == 1 && childAllele2 == 0)
			{				
				return 16; 
			}
			else if(childAllele1 == 0 && childAllele2 == 0)
			{
				return 10;
			}
			else //child is 1 / 1
			{
				return 8;
			};
		};
    };

  //Mendelian error
  return 0;//error

};

//! Returns the genotype group (1-15) a child parents trio belongs.
unsigned int CountedTrios::getGroup(const unsigned int & snpId, Subject * father, Subject * mother, Subject * child)
{

	pair<bool, bool> fatherGenotype, motherGenotype, childGenotype;
	fatherGenotype = father->getGenotype(snpId);
	motherGenotype = mother->getGenotype(snpId);
	childGenotype = child->getGenotype(snpId);

	if(motherGenotype.first == 1 && motherGenotype.second == 1)
	{
		if(fatherGenotype.first == 1 && fatherGenotype.second == 1)
		{
			if(childGenotype.first == 1 && childGenotype.second == 1)
			{
				return 1;
			};
			//else error
		}
		else if(fatherGenotype.first == 0 && fatherGenotype.second == 0)
		{
			if(childGenotype.first == 0 && childGenotype.second == 1)
			{
				return 6;
			};
			//else error
		}
		else //father is 0 / 1
		{
			if(childGenotype.first == 0 && childGenotype.second == 1)
			{
				return 3;
			}
			else if(childGenotype.first == 1 && childGenotype.second == 1)
			{
				return 2;
			};
			//else error
		};
    }
	else if(motherGenotype.first == 0 && motherGenotype.second == 0)
	{
		if(fatherGenotype.first == 1 && fatherGenotype.second == 1)
		{
			if(childGenotype.first == 0 && childGenotype.second == 1)
			{
				return 7;
			};
			//else error
		}
		else if(fatherGenotype.first == 0 && fatherGenotype.second == 0)
		{

			if(childGenotype.first == 0 && childGenotype.second == 0)
			{
				return 15;
			};
			//else error
		}
		else//father is 0 / 1
		{
			if(childGenotype.first == 0 && childGenotype.second == 1)
			{
				return 13;
			}
			else if(childGenotype.first == 0 && childGenotype.second == 0)
			{
				return 14;
			}
			//else error
		};
	}
	else //mother is 0 / 1
	{
		if(fatherGenotype.first == 1 && fatherGenotype.second == 1)
		{
			if(childGenotype.first == 0 && childGenotype.second == 1)
			{
				return 5;
			}
			else if(childGenotype.first == 1 && childGenotype.second == 1)
			{
				return 4;
			};
			//else error
		}
		else if(fatherGenotype.first == 0 && fatherGenotype.second == 0)
		{
			if(childGenotype.first == 0 && childGenotype.second == 1)
			{
				return 11;
			}
			else if(childGenotype.first == 0 && childGenotype.second == 0)
			{
				return 12;
			};
			//else error
		}
		else//father is 0 / 1
		{
			if(childGenotype.first == 0 && childGenotype.second == 1)
			{
				return 9;
			}
			else if(childGenotype.first == 0 && childGenotype.second == 0)
			{
				return 10;
			}
			else //child is 1 / 1
			{
				return 8;
			};
		};
    };

  //Mendelian error
  return 0;//error
};

//! Changes the genotype group IDs for trios by switching the allele labelling.
void CountedTrios::changeGroupIdForReverseLabels(unsigned int & groupId)
{
	if(groupId == 0) return;
	if(groupId <= 5)
	{
		groupId = 16 - groupId;
	}
	else if(groupId <= 7)
	{
		groupId = 13 - groupId;
	}
	else if(groupId <= 10)
	{
		groupId = 18 - groupId;
	}
	else if(groupId <= 15)
	{
		groupId = 16 - groupId;
	}
	else if(groupId == 16) groupId = 17; //groups 9a and 9b must switch also
	else if(groupId == 17) groupId = 16;


};

//! Adds a counted trio for phased data.
void CountedTrios::addTrioPhased(const unsigned int & snpId, const bool & reverseAlleleLabels, bool & fatherAllele1, bool & fatherAllele2, bool & motherAllele1, bool & motherAllele2, bool & childAllele1, bool & childAllele2)
{
	unsigned int groupId = getGroupPhased(snpId, fatherAllele1, fatherAllele2, motherAllele1, motherAllele2, childAllele1, childAllele2);

	if(reverseAlleleLabels) changeGroupIdForReverseLabels(groupId);
	addItemPhased(snpId, groupId);
};

//! Adds a counted trio.
void CountedTrios::addTrio(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * father, Subject * mother, Subject * child)
{
	unsigned int groupId = getGroup(snpId, father, mother, child);

	if(reverseAlleleLabels) changeGroupIdForReverseLabels(groupId);
	addItem(snpId, groupId);
};

//! Returns the genotype group (1-7) a duo belongs.
unsigned int CountedDuos::getGroup(const unsigned int & snpId, Subject * parent, Subject * child)
{

  pair<bool, bool> parentGenotype, childGenotype;
  parentGenotype = parent->getGenotype(snpId);
  childGenotype = child->getGenotype(snpId);

  if(parentGenotype.first == 1 && parentGenotype.second == 1)
    {
     
      if(childGenotype.first == 1 && childGenotype.second == 1)
        {
          return 1;
        }
      else if(childGenotype.first == 0 && childGenotype.second == 1)
	{
	  return 2;
	};
      //else error
      
    }
  else if(parentGenotype.first == 0 && parentGenotype.second == 0)
    {

      if(childGenotype.first == 0 && childGenotype.second == 0)
        {
          return 7;
        }
      else if(childGenotype.first == 0 && childGenotype.second == 1)
		{
		  return 6;
		};
      //else error

    }
	else //parent is 0 / 1
    {

      if(childGenotype.first == 0 && childGenotype.second == 0)
        {
          return 5;
        }
      else if(childGenotype.first == 0 && childGenotype.second == 1)
		{
			return 4;
		}
      else //child is 1 / 1
		{
		return 3;
		};
    };

  //Mendelian error
  //error
  return 0;
};

//! Returns the genotype group (1-7) a duo belongs, plus 4a(8) and 4b(9) for phased data.
unsigned int CountedDuos::getGroupPhased(const unsigned int & snpId, bool & parentAllele1, bool & parentAllele2, bool & childAllele1, bool & childAllele2, const bool & mother)
{

  if(parentAllele1 == 1 && parentAllele2 == 1)
    {
     
      if(childAllele1 == 1 && childAllele2 == 1)
        {
          return 1;
        }
      else if((childAllele1 == 0 && childAllele2 == 1) || (childAllele1 == 1 && childAllele2 == 0))
	{
	  return 2;
	};
      //else error
      
    }
  else if(parentAllele1 == 0 && parentAllele2 == 0)
    {

      if(childAllele1 == 0 && childAllele2 == 0)
        {
          return 7;
        }
      else if((childAllele1 == 0 && childAllele2 == 1) || (childAllele1 == 1 && childAllele2 == 0))
		{
		  return 6;
		};
      //else error

    }
  else //parent is 0 / 1
    {

		if(childAllele1 == 0 && childAllele2 == 0)
        {
			return 5;
        }
		else if(childAllele1 == 0 && childAllele2 == 1) //first allele is from the parent, as given by SHAPEIT
		{	
			if(mother) return 8; //was 4
			else return 9;
		}
		else if(childAllele1 == 1 && childAllele2 == 0) 
		{	
			if(mother) return 9; //was 4
			else return 8;
		}
		else //child is 1 / 1
		{
			return 3;
		};
    };

  //Mendelian error
  //error
  return 0;
};

//! Changes the genotype group IDs for duos by switching the allele labelling.
void CountedDuos::changeGroupIdForReverseLabels(unsigned int & groupId)
{
	if(groupId == 0) return;
	
	if(groupId < 8) groupId = 8 - groupId;	
	else if(groupId == 8) groupId = 9; //groups 4a and 4b must switch also
	else if(groupId == 9) groupId = 8;
};

//! Adds a counted duo for phased data.
void CountedDuos::addDuoPhased(const unsigned int & snpId, const bool & reverseAlleleLabels, bool & parentAllele1, bool & parentAllele2, bool & childAllele1, bool & childAllele2, const bool & mother)
{
	unsigned int groupId = getGroupPhased(snpId, parentAllele1, parentAllele2, childAllele1, childAllele2, mother);

	if(reverseAlleleLabels) changeGroupIdForReverseLabels(groupId);

	addItemPhased(snpId, groupId);
};

//! Adjusts counts so that is not inflated.
void CountedDuos::adjustPhasedDuoData(CountAlleles * countAlleles, unsigned int & noOfSNPs, const bool & mother)
{
	//should be equal for all SNPs so use the first SNP that is not zero
	map<unsigned int, map<unsigned int, unsigned int> >::iterator cdp = countedDataPhased.begin();
	unsigned int noDuos = 0;
	
	while(noDuos==0 && cdp != countedDataPhased.end())
	{
	  noDuos = cdp->second[1] + cdp->second[2] + cdp->second[3]
			  + cdp->second[5] + cdp->second[6]
			  + cdp->second[7] + cdp->second[8] + cdp->second[9];
	  cdp++;
	};

	if(noDuos == 0) return; //nothing to adjust

	//set up model for fitting
	ModelDuoAdjust * modelDuoAdj = new ModelDuoAdjust(mother);
	
	map<unsigned int, double> mafs = countAlleles->getMinorAlleleFreqs(noOfSNPs); //from all data
	bool reverseAlleles;
	bool allele1Minor;

	//loop through phased data snp by snp, use all counted data later
	for(map<unsigned int, map<unsigned int, unsigned int> >::iterator cdp = countedDataPhased.begin(); cdp != countedDataPhased.end(); ++cdp)
	{
		reverseAlleles = countAlleles->getWereAllelesReversed(cdp->first);
		allele1Minor = countAlleles->isSnpAllele1Minor(cdp->first);

		//add minor alleles from father and mother
		if((allele1Minor && !reverseAlleles) || (!allele1Minor && reverseAlleles))
		{
			modelDuoAdj->addFromFatherData(cdp->first, cdp->second[9]);
			modelDuoAdj->addFromMotherData(cdp->first, cdp->second[8]);
		}
		else
		{
			modelDuoAdj->addFromFatherData(cdp->first, cdp->second[8]);
			modelDuoAdj->addFromMotherData(cdp->first, cdp->second[9]);		
		};
	};
	
	modelDuoAdj->setNoOfDuos(noDuos);
	modelDuoAdj->addMAFData(mafs);

	//set up initial parameters
	modelDuoAdj->setParameter(1, 0);       
	modelDuoAdj->setParameter(2, 0);
	modelDuoAdj->setParameter(3, 0);
	modelDuoAdj->setParameter(4, 0);

	FindFit findFit(modelDuoAdj);

	//find best adjustment
	double eval;

	double accuracy = 1e-6;
	unsigned int maxIterations = 50;
	bool fitOK = findFit.newtonsMethod(eval, accuracy, maxIterations);

	if(!fitOK) //reset parameters
	{
		modelDuoAdj->setParameter(1, 0);       
		modelDuoAdj->setParameter(2, 0);
		modelDuoAdj->setParameter(3, 0);
		modelDuoAdj->setParameter(4, 0);

		out("Warning! Failed to adjust ");
		if(mother) out("case/mother duo cell counts!\n\n"); else out("case/father duo cell counts!\n\n");
	};	

	//now adjust the counts
	double beta0 = modelDuoAdj->getParameter(1);
	double beta1 = modelDuoAdj->getParameter(2);
	double beta2 = modelDuoAdj->getParameter(3);
	double beta3 = modelDuoAdj->getParameter(4);

	if(mother) out("Case/mother duo adjustment parameters:\n"); else out("Case/father duo adjustment parameters:\n");
	out("beta0 = "); out(beta0); out("\n");
	out("beta1 = "); out(beta1); out("\n");
	out("beta2 = "); out(beta2); out("\n");
	out("beta3 = "); out(beta3); out("\n\n");

	double adjustment, mafSq;
	map<unsigned int, double> countedDataAdj; //group ID, adj count
	
	map<unsigned int, double>::const_iterator mafi = mafs.begin();

	for(map<unsigned int, map<unsigned int, unsigned int> >::iterator cdpa = countedDataPhased.begin(); cdpa != countedDataPhased.end(); ++cdpa)
	{
		//ensure correct minor allele freq is used, skip unused SNPs
		while(mafi->first != cdpa->first && mafi != mafs.end()) mafi++;
	
		//should not need this bit
		if(mafi->first != cdpa->first)
		{
			mafi = mafs.find(cdpa->first);
			if(mafi == mafs.end())
			{
				outErr("Error: SNP minor allele value missing for SNP "); outErr(cdpa->first); outErr(" unable to adjust phased duos!\n\n");
				exit(1);
			};
		};

	
		//set adjustment value		
		mafSq = (mafi->second)*(mafi->second);
		adjustment = (double)(noDuos)*(beta0 + beta1*(mafi->second) + beta2*mafSq + beta3*mafSq*(mafi->second));

		reverseAlleles = countAlleles->getWereAllelesReversed(cdpa->first); //that is, risk allele is the major allele if true
		allele1Minor = countAlleles->isSnpAllele1Minor(cdpa->first);

		//make sure adjustment subtracts from cell with number of minor alleles inherited from father, *not* risk allele
		if((allele1Minor && !reverseAlleles) || (!allele1Minor && reverseAlleles)) adjustment = -adjustment;

		//set new adjusted estimated values
		if(cdpa->second[8] - adjustment < 0)
		{
			countedDataAdj[8] = 0;
			countedDataAdj[9] = cdpa->second[8] + cdpa->second[9];
		}
		else if(cdpa->second[9] + adjustment < 0)
		{
			countedDataAdj[8] = cdpa->second[8] + cdpa->second[9];
			countedDataAdj[9] = 0;
		}
		else
		{
			countedDataAdj[8] = cdpa->second[8] - adjustment; //from father, update for each SNP
			countedDataAdj[9] = cdpa->second[9] + adjustment; //from mother
		};

		countedDataPhasedAdjusted[cdpa->first] = countedDataAdj; //add adjusted data for each SNP
		++mafi;
	};

	delete modelDuoAdj;
};

//! Adds a duo.
void CountedDuos::addDuo(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * parent, Subject * child)
{
  unsigned int groupId = getGroup(snpId, parent, child);

  if(reverseAlleleLabels) changeGroupIdForReverseLabels(groupId);
  addItem(snpId, groupId);

};

//! Returns the genotype group (1-9) a pair of parents belongs.
unsigned int CountedParents::getGroup(const unsigned int & snpId, Subject * father, Subject * mother)
{
  pair<bool, bool> fatherGenotype, motherGenotype;
  fatherGenotype = father->getGenotype(snpId);
  motherGenotype = mother->getGenotype(snpId);

  if(motherGenotype.first == 1 && motherGenotype.second == 1)
    {
      if(fatherGenotype.first == 1 && fatherGenotype.second == 1)
	{
	  return 1;
	}
      else if(fatherGenotype.first == 0 && fatherGenotype.second == 0)
	{
	  return 3;
	}
      else //father is 0 / 1
	{
	  return 2;
	};

    }
  else if(motherGenotype.first == 0 && motherGenotype.second == 0)
    {

      if(fatherGenotype.first == 1 && fatherGenotype.second == 1)
	{
	  return 7;
	}
      else if(fatherGenotype.first == 0 && fatherGenotype.second == 0)
	{
	  return 9;
	}
      else//father is 0 / 1
	{
	  return 8;
	};

    }
  else //mother is 0 / 1
    {
      if(fatherGenotype.first == 1 && fatherGenotype.second == 1)
	{
	  return 4;
	}
      else if(fatherGenotype.first == 0 && fatherGenotype.second == 0)
	{
	  return 6;
	}
      else//father is 0 / 1
	{
	  return 5;
	};

    };

  return 0;//error
};

//! Changes the genotype group IDs for parents by switching the allele labelling.
void CountedParents::changeGroupIdForReverseLabels(unsigned int & groupId)
{
	if(groupId == 0) return;
	
	groupId = 10 - groupId;
	
};

//! Adds a pair of parents.
void CountedParents::addParents(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * father, Subject * mother)
{

  unsigned int groupId = getGroup(snpId, father, mother);

  if(reverseAlleleLabels) changeGroupIdForReverseLabels(groupId);
  addItem(snpId, groupId);
};

//! Returns the genotype group (1-3) a subject belongs.
unsigned int CountedSingleSubject::getGroup(const unsigned int & snpId, Subject * subject)
{
	pair<bool, bool> subjectGenotype;
	subjectGenotype = subject->getGenotype(snpId);

	if(subjectGenotype.first == 1 && subjectGenotype.second == 1)
	{
		return 1;
	}
	else if(subjectGenotype.first == 0 && subjectGenotype.second == 0)
	{
		return 3;
	}
	else//subject is 0 / 1
	{
		return 2;
	};

	return 0;
};

//! Adds a subject.
void CountedSingleSubject::addSubject(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * subject)
{

  unsigned int groupId = getGroup(snpId, subject);

  if(reverseAlleleLabels) changeGroupIdForReverseLabels(groupId);
  addItem(snpId, groupId);
};

//! Changes the genotype group IDs for parents by switching the allele labelling.
void CountedSingleSubject::changeGroupIdForReverseLabels(unsigned int & groupId)
{
	if(groupId == 0) return;
	
	groupId = 4 - groupId;
	
};

//! Returns a subject in a pedigree
Subject * Pedigree::getSubject(const unsigned int & id) const
{
	map<unsigned int, Subject *>::const_iterator s = subjects.find(id);
	if(s != subjects.end()) return s->second;

	outErr("Error: subject not found in pedigree!\n");
	exit(1);	       
};

//! Returns a subject in a pedigree and returns if it exists. 
Subject * Pedigree::getSubject(const unsigned int & id, bool & exists) const
{
	map<unsigned int, Subject *>::const_iterator s = subjects.find(id);
	if(s != subjects.end())
	{
		exists = true;
		return s->second;
	};

	exists = false;
	return 0;       
};

//! Checks if a subject is a child of a case parent trio in a pedigree.
void Pedigree::findCaseParentTrio(Trio & aTrio, const unsigned int & snpId, bool & found, Subject * subject)
{
	bool fatherExists = false;
	bool motherExists = false;
	
	//check if subject is affected, has genotype data and has not been removed if allowing extra trios 
	if(subject->getAffected() && subject->genotypeDataPresent(snpId) && subject->isNotRemoved())
	{
		//check if parents exist and have genotype data
		if(subject->getFatherId() != 0 && subject->getMotherId() != 0)
		{
			Subject * father = getSubject(subject->getFatherId(), fatherExists);
			Subject * mother = getSubject(subject->getMotherId(), motherExists);
			if(fatherExists && motherExists && father->isNotRemoved() && mother->isNotRemoved())
			{
				if(father->genotypeDataPresent(snpId) && mother->genotypeDataPresent(snpId))
				{ 
					aTrio.father = father;
					aTrio.mother = mother;
					aTrio.child = subject;
					found = true;
				};
			};
		};
	};

};

//! Loops thro' subjects in pedigree to find a case parent trio.
Trio Pedigree::findCaseParentTrio(const unsigned int & snpId, bool & found)
{
	found = false;
	Trio aTrio;	

	//loop through subjects looking for a case which has parents
	for(map<unsigned int, Subject *>::const_iterator s = subjects.begin(); (s != subjects.end() && !found); ++s)
	{
		 findCaseParentTrio(aTrio, snpId, found, s->second);	
	};

	return aTrio;
};

//! Adds count to suitable trios for a pedigree
void addToChildrenCountForTrioOrDuo(map<unsigned int, map<unsigned int, unsigned int> > & childrenCount, const unsigned int & pedId, const unsigned int & childId)
{
	map<unsigned int, map<unsigned int, unsigned int> >::iterator cctd = childrenCount.find(pedId);
	if(cctd != childrenCount.end())
	{
		map<unsigned int, unsigned int>::iterator ctd = cctd->second.find(childId);
		if(ctd != cctd->second.end()) ctd->second++;
		else ctd->second = 1; //this shouldn't be possible
	}
	else
	{
		map<unsigned int, unsigned int> newCount;
		newCount[childId] = 1;
		childrenCount[pedId] = newCount;
	};
};

//! Loop thro' subjects and keep track of all valid case/parent trios
bool Pedigree::findPooHapCaseParentTrios(const unsigned int & snpId, map<unsigned int, map<unsigned int, unsigned int> > & childrenCountForTrio)
{
	bool found = false;
	bool anyFound = false;
	Trio aTrio;	

	//loop through subjects looking for which cases have parents
	for(map<unsigned int, Subject *>::const_iterator s = subjects.begin(); s != subjects.end(); ++s)
	{
		found = false;
		findCaseParentTrio(aTrio, snpId, found, s->second);	//is this subject a child of a case/parent trio?

		if(found) {addToChildrenCountForTrioOrDuo(childrenCountForTrio, pedId, s->second->getId()); anyFound = true;};
	};

	return anyFound;
};

//! Loop thro' subjects and keep track of all valid case/mother duos
bool Pedigree::findPooHapCaseMotherDuos(const unsigned int & snpId, map<unsigned int, map<unsigned int, unsigned int> > & childrenCountForMotherDuo)
{
	bool found = false;
	bool anyFound = false;
	pair<Subject *, Subject *> aCaseMotherDuo;

	//loop through subjects looking for which cases have parents
	for(map<unsigned int, Subject *>::const_iterator s = subjects.begin(); s != subjects.end(); ++s)
	{
		found = false;		
		findCaseMotherDuo(aCaseMotherDuo, snpId, found, s->second);	//is this subject a child of a case/mother duo?

		if(found) {addToChildrenCountForTrioOrDuo(childrenCountForMotherDuo, pedId, s->second->getId()); anyFound = true;};
	};

	return anyFound;
};

//! Loop thro' subjects and keep track of all valid case/father duos
bool Pedigree::findPooHapCaseFatherDuos(const unsigned int & snpId, map<unsigned int, map<unsigned int, unsigned int> > & childrenCountForFatherDuo)
{
	bool found = false;
	bool anyFound = false;
	pair<Subject *, Subject *> aCaseFatherDuo;

	//loop through subjects looking for which cases have parents
	for(map<unsigned int, Subject *>::const_iterator s = subjects.begin(); s != subjects.end(); ++s)
	{
		found = false;		
		findCaseFatherDuo(aCaseFatherDuo, snpId, found, s->second);	//is this subject a child of a case/father duo?

		if(found) {addToChildrenCountForTrioOrDuo(childrenCountForFatherDuo, pedId, s->second->getId()); anyFound = true;};
	};

	return anyFound;
};

//! Checks if a subject is a child of a case mother duo in a pedigree.
void Pedigree::findCaseMotherDuo(pair<Subject *, Subject *> & aCaseMotherDuo, const unsigned int & snpId, bool & found, Subject * subject)
{
	bool motherExists = false;

	//check if we have a case with genotype data
	if(subject->getAffected() && subject->genotypeDataPresent(snpId))
	{
		//check if case has mother with genotype data
		if(subject->getMotherId() != 0)
		{
			Subject * mother = getSubject(subject->getMotherId(), motherExists);
			if(motherExists)
			{
				if(mother->genotypeDataPresent(snpId))
				{
					aCaseMotherDuo = make_pair(subject, mother);
					found = true;
				};
			};
		};
	};
	
};

//! Loops thro' subjects in pedigree to find a case mother duo.
pair<Subject *, Subject *> Pedigree::findCaseMotherDuo(const unsigned int & snpId, bool & found)
{
	found = false;
	pair<Subject *, Subject *> aCaseMotherDuo;

	//loop through subjects looking for a case with a mother
	for(map<unsigned int, Subject *>::const_iterator s = subjects.begin(); (s != subjects.end() && !found); ++s)
	{
		findCaseMotherDuo(aCaseMotherDuo, snpId, found, s->second);
	};

	return aCaseMotherDuo;
};

//! Checks if a subject is a child of a case father duo in a pedigree.
void Pedigree::findCaseFatherDuo(pair<Subject *, Subject *> & aCaseFatherDuo, const unsigned int & snpId, bool & found, Subject * subject)
{
	bool fatherExists = false;
	
	//check if we have a case with genotype data
	if(subject->getAffected() && subject->genotypeDataPresent(snpId))
	{
		//check if case has mother with genotype data
		if(subject->getFatherId() != 0)
		{
			Subject * father = getSubject(subject->getFatherId(), fatherExists);
			if(fatherExists)
			{
				if(father->genotypeDataPresent(snpId))
				{
					aCaseFatherDuo = make_pair(subject, father);
					found = true;
				};
			};
		};
	};
	
};

//! Loops thro' subjects in pedigree to find a case father duo.
pair<Subject *, Subject *> Pedigree::findCaseFatherDuo(const unsigned int & snpId, bool & found)
{
	found = false;
	pair<Subject *, Subject *> aCaseFatherDuo;

	//loop through subjects looking for a case with a father
	for(map<unsigned int, Subject *>::const_iterator s = subjects.begin(); (s != subjects.end() && !found); ++s)
	{
		findCaseFatherDuo(aCaseFatherDuo, snpId, found, s->second);
	};

	return aCaseFatherDuo;
};

//! Checks if a subject is a case in a pedigree.
void Pedigree::findCase(Subject *& aCase, const unsigned int & snpId, bool & found, Subject * subject)
{
	if(subject->getAffected() && subject->genotypeDataPresent(snpId))
	{
		aCase = subject;
		found = true; 
	};
	
};

//! Loops thro' subjects in pedigree to find a case.
Subject * Pedigree::findCase(const unsigned int & snpId, bool & found)
{
	found = false;
	Subject * aCase = 0;

	//loop through subjects looking for a case with genotype data
	for(map<unsigned int, Subject *>::const_iterator s = subjects.begin(); (s != subjects.end() && !found); ++s)
	{
		findCase(aCase, snpId, found, s->second);
	};

	return aCase;
};

//! Checks if a subject is a case with parents in a pedigree.
void Pedigree::findParentsOfCase(pair<Subject *, Subject *> & parents, const unsigned int & snpId, bool & found, Subject * subject)
{
	bool fatherExists = false;
	bool motherExists = false;
	
	//check if subject is affected, but do not check genotype data
	if(subject->getAffected())
	{
		//check if parents exist and have genotype data
		if(subject->getFatherId() != 0 && subject->getMotherId() != 0)
		{
			Subject * father = getSubject(subject->getFatherId(), fatherExists);
			Subject * mother = getSubject(subject->getMotherId(), motherExists);
			if(fatherExists && motherExists)
			{
				if(father->genotypeDataPresent(snpId) && mother->genotypeDataPresent(snpId))
				{ 
					parents = make_pair(father, mother);
					found = true;
				};
			};
		};
	};
	
};

//! Loops thro' subjects in pedigree to find parents of a case.
pair<Subject *, Subject *> Pedigree::findParentsOfCase(const unsigned int & snpId, bool & found)
{
	found = false;
	pair<Subject *, Subject *> parents;//father, mother

	//loop through subjects looking for a case which has parents
	for(map<unsigned int, Subject *>::const_iterator s = subjects.begin(); (s != subjects.end() && !found); ++s)
	{
		findParentsOfCase(parents, snpId, found, s->second);
	};

	return parents;
};

//! Checks if a subject is a case with a mother in a pedigree.
void Pedigree::findMotherOfCase(Subject *& mother, const unsigned int & snpId, bool & found, Subject * subject)
{	
	bool motherExists = false;
		
	//check if subject is affected, but do not check genotype data
	if(subject->getAffected())
	{
		//check if mother exists and has genotype data
		if(subject->getMotherId() != 0)
		{
			mother = getSubject(subject->getMotherId(), motherExists);
			if(motherExists)
			{
				if(mother->genotypeDataPresent(snpId))
				{ 
					found = true;
				};
			};
		};
	};
	
};

//! Loops thro' subjects in pedigree to find a mother of a case.
Subject * Pedigree::findMotherOfCase(const unsigned int & snpId, bool & found)
{
	found = false;
	Subject * mother = 0;

	//loop through subjects looking for a case with a mother
	for(map<unsigned int, Subject *>::const_iterator s = subjects.begin(); (s != subjects.end() && !found); ++s)
	{
		findMotherOfCase(mother, snpId, found, s->second);
	};

	return mother;
};

//! Checks if a subject is a case with a father in a pedigree.
void Pedigree::findFatherOfCase(Subject *& father, const unsigned int & snpId, bool & found, Subject * subject)
{
	bool fatherExists = false;
	
	//check if subject is affected, but do not check genotype data
	if(subject->getAffected())
	{
		//check if father exists and has genotype data
		if(subject->getFatherId() != 0)
		{
			father = getSubject(subject->getFatherId(), fatherExists);
			if(fatherExists)
			{
				if(father->genotypeDataPresent(snpId))
				{ 
					found = true;
				};
			};
		};
	};
	
};

//! Loops thro' subjects in pedigree to find a father of a case.
Subject * Pedigree::findFatherOfCase(const unsigned int & snpId, bool & found)
{
	found = false;
	Subject * father = 0;

	//loop through subjects looking for a case with a father
	for(map<unsigned int, Subject *>::const_iterator s = subjects.begin(); (s != subjects.end() && !found); ++s)
	{
		findFatherOfCase(father, snpId, found, s->second);	
	};

	return father;
};

//! Checks if a subject is a control with parents in a pedigree.
void Pedigree::findParentsOfControl(Trio & aTrio, const unsigned int & snpId, bool & found, Subject * subject)
{
	bool fatherExists = false;
	bool motherExists = false;
	
	//check if subject is not removed if allowing for extra parents of controls
	if(subject->isNotRemoved())
	{
		//check if parents exist and have genotype data
		if(subject->getFatherId() != 0 && subject->getMotherId() != 0)
		{
			Subject * father = getSubject(subject->getFatherId(), fatherExists);
			Subject * mother = getSubject(subject->getMotherId(), motherExists);
			if(fatherExists && motherExists && father->isNotRemoved() && mother->isNotRemoved())
			{
				if(father->genotypeDataPresent(snpId) && mother->genotypeDataPresent(snpId))
				{ 						
					aTrio.father = father;
					aTrio.mother = mother;
					aTrio.child = subject;
					found = true;
				};
			};
		};
	};
	
};

//! Loops thro' subjects in pedigree to find parents of a control.
Trio Pedigree::findParentsOfControl(const unsigned int & snpId, bool & found)
{
	found = false;
	Trio aTrio;

	//loop through subjects looking for a subject with unknown genotype with parents
	for(map<unsigned int, Subject *>::const_iterator s = subjects.begin(); (s != subjects.end() && !found); ++s)
	{
		findParentsOfControl(aTrio, snpId, found, s->second);
	};

	return aTrio;
};

//! Checks if a subject is a control a mother in a pedigree.
void Pedigree::findControlMotherDuo(pair<Subject *, Subject *> & aControlMotherDuo, const unsigned int & snpId, bool & found, Subject * subject)
{
	bool motherExists = false;

	//check if we have a control with genotype data
	if(!(subject->getAffected()) && subject->genotypeDataPresent(snpId))
	{
		//check if control has mother with genotype data
		if(subject->getMotherId() != 0)
		{
			Subject * mother = getSubject(subject->getMotherId(), motherExists);

			if(motherExists)
			{
				if(mother->genotypeDataPresent(snpId))
				{
					aControlMotherDuo = make_pair(subject, mother);
					found = true;
				};
			};
		};
	};
	
};

//! Loops thro' subjects in pedigree to find a control mother duo.
pair<Subject *, Subject *> Pedigree::findControlMotherDuo(const unsigned int & snpId, bool & found)
{
	found = false;
	pair<Subject *, Subject *> aControlMotherDuo;

	//loop through subjects looking for a control with a mother
	for(map<unsigned int, Subject *>::const_iterator s = subjects.begin(); (s != subjects.end() && !found); ++s)
	{
		findControlMotherDuo(aControlMotherDuo, snpId, found, s->second);	
	};

	return aControlMotherDuo;
};

//! Checks if a subject is a control a father in a pedigree.
void Pedigree::findControlFatherDuo(pair<Subject *, Subject *> & aControlFatherDuo, const unsigned int & snpId, bool & found, Subject * subject)
{
	bool fatherExists = false;
	
	//check if we have a control with genotype data
	if(!(subject->getAffected()) && subject->genotypeDataPresent(snpId))
	{
		//check if control has father with genotype data
		if(subject->getFatherId() != 0)
		{
			Subject * father = getSubject(subject->getFatherId(), fatherExists);
			if(fatherExists)
			{
				if(father->genotypeDataPresent(snpId))
				{
					aControlFatherDuo = make_pair(subject, father);
					found = true;
				};
			};
		};
	};
	
};

//! Loops thro' subjects in pedigree to find a control father duo.
pair<Subject *, Subject *> Pedigree::findControlFatherDuo(const unsigned int & snpId, bool & found)
{
	found = false;
	pair<Subject *, Subject *> aControlFatherDuo;

	//loop through subjects looking for a control with a father
	for(map<unsigned int, Subject *>::const_iterator s = subjects.begin(); (s != subjects.end() && !found); ++s)
	{
		findControlFatherDuo(aControlFatherDuo, snpId, found, s->second);
	};

	return aControlFatherDuo;
};

//! Checks if a subject is a control in a pedigree.
void Pedigree::findControl(Subject *& control, const unsigned int & snpId, bool & found, Subject * subject)
{
	
	//check if we have a control with genotype data
	if(!(subject->getAffected()) && subject->genotypeDataPresent(snpId))
	{		
		control = subject;
		found = true;		
	};
	
};

//! Loops thro' subjects in pedigree to find a control.
Subject * Pedigree::findControl(const unsigned int & snpId, bool & found)
{
	found = false;
	Subject * control = 0;

	//loop through subjects looking for a control
	for(map<unsigned int, Subject *>::const_iterator s = subjects.begin(); (s != subjects.end() && !found); ++s)
	{
		findControl(control, snpId, found, s->second);
	};
	
	return control;
};

//! Unmark all subjects in pedigree for when adding extra trios.
void Pedigree::restoreAllSubjectsToPedigree()
{
	for(map<unsigned int, Subject *>::const_iterator s = subjects.begin(); s != subjects.end(); ++s)
	{
		s->second->restoreToPedigree();
	};
};

