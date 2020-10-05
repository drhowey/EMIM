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


/*! \file ProcessData.h
    \brief This file contains classes and methods for processing the pedigree data into separate genotype groups.
    
  
*/

#ifndef __PROCESSDATA
#define __PROCESSDATA

#include <string>
#include <map>
#include <set>
#include <list>
#include <iostream>
#include <fstream>

#include "Data.h"
#include "main.h"

#ifdef USING_GZIP
#include <gzstream.h>
#endif


 //! namespace std initiates the "std" or "standard" namespace.
using namespace std;

//! Contains a summary of the counted genotypes groups from pedigrees.
class AllCountedData
{
private:

	CountedTrios caseParentTrios;
	CountedDuos caseFatherDuos;
	CountedDuos caseMotherDuos;
	CountedDuos controlFatherDuos;
	CountedDuos controlMotherDuos;
	CountedParents parentsOfCases;
	CountedParents parentsOfControls;
	CountedSingleSubject fathersOfCases;
	CountedSingleSubject mothersOfCases;
	CountedSingleSubject cases;
	CountedSingleSubject controls;
	
public:

	//! Create an object to hold all of the processed pedigree infomation.
	/*!
	    The output file name for each pedigree subgroup type is set and
	 the number of different genotype groups for each. 
	 */
	AllCountedData(const string & dir, const string & en, const unsigned int & sso, bool & ph) : caseParentTrios("caseparenttrios", 15, "case parent trios", dir, en, sso, ph), caseFatherDuos("casefatherduos", 7, "case father duos", dir, en, sso, ph),  caseMotherDuos("casemotherduos", 7, "case mother duos", dir, en, sso, ph),
		controlFatherDuos("confatherduos", 7, "control father duos", dir, en, sso), controlMotherDuos("conmotherduos", 7, "control mother duos", dir, en, sso), parentsOfCases("caseparents", 9, "case parents", dir, en, sso),
		parentsOfControls("conparents", 9, "control parents", dir, en, sso), fathersOfCases("casefathers", 3, "case fathers", dir, en, sso), mothersOfCases("casemothers", 3, "case mothers", dir, en, sso), cases("cases", 3, "cases", dir, en, sso),
		controls("cons", 3, "controls", dir, en, sso) {};

	void addCaseParentTrioPhased(const unsigned int & snpId, const bool & reverseAlleleLabels, bool & fatherAllele1, bool & fatherAllele2, bool & motherAllele1, bool & motherAllele2, bool & childAllele1, bool & childAllele2);
	void addCaseMotherDuoPhased(const unsigned int & snpId, const bool & reverseAlleleLabels, bool & motherAllele1, bool & motherAllele2, bool & childAllele1, bool & childAllele2);
	void addCaseFatherDuoPhased(const unsigned int & snpId, const bool & reverseAlleleLabels, bool & motherAllele1, bool & motherAllele2, bool & childAllele1, bool & childAllele2);
	void addCaseParentTrio(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * father, Subject * mother, Subject * child);
	void addCaseFatherDuo(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * father, Subject * child);
	void addCaseMotherDuo(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * mother, Subject * child);
	void addControlFatherDuo(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * father, Subject * child);
	void addControlMotherDuo(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * mother, Subject * child);
	void addParentsOfCase(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * father, Subject * mother);
	void addParentsOfControl(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * father, Subject * mother);
	void addFatherOfCase(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * father);
	void addMotherOfCase(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * mother);
	void addCase(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * subject);
	void addControl(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * subject);
	void outputResults(const unsigned int & noOfSnps, const unsigned int & totalNoOfPedigrees, const unsigned int & splitSNPOutput);
	void outputFileSummaries(const unsigned int & noOfSnps, const unsigned int & totalNoOfPedigrees);
	void outputResultOneSnp(const unsigned int & snpId, const unsigned int & splitSNPOutput, bool & pooHaps);
	void outputParameterFile(ofstream & fileOutPut, const unsigned int & noOfSnps, const bool & useCPG);
	void setTotalSamples(unsigned int & totSams) {caseParentTrios.setTotalSamples(totSams); caseMotherDuos.setTotalSamples(totSams); caseFatherDuos.setTotalSamples(totSams);};
	void setNoOfTrioGroups(const unsigned int & ng) {caseParentTrios.setNoOfGroups(ng);};
	void setNoOfDuoGroups(const unsigned int & ng) {caseMotherDuos.setNoOfGroups(ng); caseFatherDuos.setNoOfGroups(ng);};
	void adjustPhasedDuoData(CountAlleles * countAlleles, unsigned int & noOfSNPs) {caseMotherDuos.adjustPhasedDuoData(countAlleles, noOfSNPs, true); caseFatherDuos.adjustPhasedDuoData(countAlleles, noOfSNPs, false);};
	void setAdjustPhasedData(bool & adj) {caseMotherDuos.setAdjustedCounts(adj); caseFatherDuos.setAdjustedCounts(adj);};
};

//! Contains a list of all pedigrees.
class AllPedigreeData
{
private:
	

public:

	map<unsigned int, Pedigree *> pedigrees; //!< pedigree ID, pedigree

	AllPedigreeData() : pedigrees() {};

	//! Delete all the pedigrees, these belong to this class
	~AllPedigreeData()
	{
		for(map<unsigned int, Pedigree *>::iterator p = pedigrees.begin(); p != pedigrees.end(); ++p)
		{
			delete p->second;
		};
	};

	void addSubjectToPedigree(unsigned int & subjectId, Subject * subject, unsigned int & pedigreeId);
	void outputSummary();
	void restoreAllSubjectsToAllPedigrees();
	void subjectExistsAndAffected(const unsigned int & pedId, const unsigned int & subjectId, MapIds & subjectIds);  
};

//! Used to organise the processing of files and data.
class ProcessData
{
private:
	bool extraAffectedTrios, extraUnaffectedTrios, childGenotype, childTrend, motherGenotype,
		motherTrend, imprintingMaternal, imprintingPaternal,
		imprintingMaternalWeinberg, imprintingPaternalWeinberg, estimateAlleleFreq, useMajorAlleleAsRisk;
	unsigned int splitSNPOutput;
	string outputDirectory;
	string endName;
	string probandFileName;
	string riskAlleleInFileName;
	string riskAlleleOutFileName;
	AllPedigreeData allPedigreeData;
	map<unsigned int, set<unsigned int> > probandSubjectIds; //pedigree id, subject ids
	AllCountedData allCountedData;
	map<unsigned int, string> snpNames;
	map<unsigned int, string> riskAlleleNames;//snp No., risk allele name
	unsigned int noOfSnps;
	unsigned int pedigreesNotCounted;
	string pedigreeName, subjectIdName, fatherIdName, motherIdName, sexIdName, affectedIdName;
	string allele1Name, allele2Name;
	unsigned int fileType; // 1 = .ped, 2 = .bed, 3 = .gzip
	MapSnpAlleleIds snpAlleleIds;
	bool allele1Id, allele2Id;	
	ifstream readPedigree;
	ifstream readBinaryGenotypeData;
	string bimFileName;

#ifdef USING_GZIP
	igzstream readGzipPedigree;
#endif
	ofstream markersFile;
	bool snpMajor;
	unsigned int one;
	unsigned int aBit;
	unsigned int bitCount;
	char buffer[1];
	unsigned int allele1, allele2;
	list<Subject *> orderListOfSubjects; //to use for binary SNP data

	//parent of origin using haplotye estimation
	bool pooHaps;
	bool adjustedCounts;
	list<unsigned int> includedSNPsInShapeIT;
	string hapMCMCOptions;
	string hapModelOptions;
	string otherShapeitOtherOptions;
	string threadOptions;
	string shapeitCommand;
	unsigned int noSamples;
	double missThres;

public:
	ProcessData(bool xa, bool xu, bool cg, bool ct, bool mg, bool mt, bool im, bool ip, bool imw,
		bool ipw, bool eaf, bool uma, unsigned int so, string od, string en, string pb, string rain, string raout, bool ph, bool adj, string & sic, string & to, string & hmcmco, string & hmodo, string & hotho, unsigned int & ns, double & mth) : extraAffectedTrios(xa),
		extraUnaffectedTrios(xu), childGenotype(cg), childTrend(ct), motherGenotype(mg),
		motherTrend(mt), imprintingMaternal(im), imprintingPaternal(ip),
		imprintingMaternalWeinberg(imw), imprintingPaternalWeinberg(ipw),
		estimateAlleleFreq(eaf), useMajorAlleleAsRisk(uma), splitSNPOutput(so), outputDirectory(od), endName(en),
		probandFileName(pb), riskAlleleInFileName(rain), riskAlleleOutFileName(raout), pooHaps(ph), adjustedCounts(adj), allPedigreeData(), probandSubjectIds(), allCountedData(od, en, so, ph),
		shapeitCommand(sic), threadOptions(to), hapMCMCOptions(hmcmco), hapModelOptions(hmodo), otherShapeitOtherOptions(hotho), noSamples(ns), missThres(mth)
	{one = '\1'; if(pooHaps) {allCountedData.setNoOfTrioGroups(17); allCountedData.setNoOfDuoGroups(9); allCountedData.setAdjustPhasedData(adj);}; };

	~ProcessData() {};

	void process(string & fileName, string & mapFileName);
	void createMarkerFiles(const unsigned int & splitSNPOutput, string & markersFile, string & outputDirectory, string & endName);
	void createResultsFiles(string & outputDirectory, string & endName);
	void outputOptionInfo(string & fileName, string & mapFileName);
	void setNoOfSnps(string & fileName, string & mapFileName);
	void analyseData(string & fileName);
	bool getReverseAlleles(CountAlleles *& countAlleles, const unsigned int & snpNo) const;
	void setRiskAlleles();
	bool addProcessedPedigree(const unsigned int & snpId, const bool & reverseAlleleLabels, Pedigree * pedigree, AllCountedData & allCountedData,
		map<unsigned int, set<unsigned int> > & probandSubjectIds);
	void processPooHaps(string & fileName, CountAlleles * countAlleles);
	bool findPooHapsTrio(const unsigned int & snpId, map<unsigned int, map<unsigned int, unsigned int> > & childrenCountForTrio, const bool & reverseAlleleLabels,
		Pedigree * pedigree, map<unsigned int, set<unsigned int> > & probandSubjectIds);
	bool findPooHapsDuo(const unsigned int & snpId, map<unsigned int, map<unsigned int, unsigned int> > & childrenCountForMotherDuo, const bool & reverseAlleleLabels,
		Pedigree * pedigree, map<unsigned int, set<unsigned int> > & probandSubjectIds, const bool & motherDuo);
	void outputPooHapTrioFamFile(const unsigned int & pedId, Pedigree * pedigree, const unsigned int & childId, ofstream & outFamPooFile);
	void outputPooHapMotherDuoFamFile(const unsigned int & pedId, Pedigree * pedigree, const unsigned int & childId, ofstream & outFamPooFile);
	void outputPooHapFatherDuoFamFile(const unsigned int & pedId, Pedigree * pedigree, const unsigned int & childId, ofstream & outFamPooFile);
	void outputPooHapBedFile(string & fileName, CountAlleles * countAlleles, map<unsigned int, unsigned int> & pooHapTrios, map<unsigned int, unsigned int> & pooHapMotherDuos, map<unsigned int, unsigned int> & pooHapFatherDuos, AllPedigreeData & allPedigreeData, set<unsigned int> & includedSNPs, ofstream & outBedPooTrioFile);
	void outputPooHapMapOrBimFile(string & fileName, string & inputFilename, set<unsigned int> & includedSNPs);
	bool addProcessedPedigreeProbandSubject(bool & found, bool & foundTrio, bool & foundParentsCon, unsigned int & probandSubjectId,
		const unsigned int & snpId, const bool & reverseAlleleLabels, Pedigree * pedigree, AllCountedData & allCountedData);
	void getSubjectNamesFromFile();
	void getAlleleNamesFromFile();
	void setEstimateAlleleFreq(bool & eaf) {estimateAlleleFreq = eaf;};
	void addGenotypeDataToSubject(Subject * subject, const unsigned int & snpId, CountAlleles * countAlleles);
	void addGenotypeDataToSubjectPed(Subject * subject, const unsigned int & snpId, CountAlleles * countAlleles);
	void addGenotypeDataToSubjectBinary(Subject * subject, const unsigned int & snpId, CountAlleles * countAlleles);
	void outputParameterFile(AllCountedData & allCountedData, const unsigned int & noOfSnps) const;
	void outputRiskAlleles(CountAlleles *& countAlleles) const;

};


#endif

