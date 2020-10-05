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


/*! \file ProcessData.cpp
    \brief This file contains classes and methods for processing the pedigree data into separate genotype groups.
    
  
*/

#include <string>
#include <map>
#include <set>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <math.h>

#include "ProcessData.h"
#include "main.h"

#ifdef USING_GZIP
#include <gzstream.h>
#endif

using namespace std; // initiates the "std" or "standard" namespace

//! Adds a case parent trio for phased data (using SHAPEIT).
void AllCountedData::addCaseParentTrioPhased(const unsigned int & snpId, const bool & reverseAlleleLabels, bool & fatherAllele1, bool & fatherAllele2, bool & motherAllele1, bool & motherAllele2, bool & childAllele1, bool & childAllele2)
{
	caseParentTrios.addTrioPhased(snpId, reverseAlleleLabels, fatherAllele1, fatherAllele2, motherAllele1, motherAllele2, childAllele1, childAllele2);
};

//! Adds a case parent trio.
void AllCountedData::addCaseParentTrio(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * father, Subject * mother, Subject * child)
{
  caseParentTrios.addTrio(snpId, reverseAlleleLabels, father, mother, child);
};

//! Adds a case/father duo for phased data (using SHAPEIT).
void AllCountedData::addCaseFatherDuoPhased(const unsigned int & snpId, const bool & reverseAlleleLabels, bool & fatherAllele1, bool & fatherAllele2, bool & childAllele1, bool & childAllele2)
{
	caseFatherDuos.addDuoPhased(snpId, reverseAlleleLabels, fatherAllele1, fatherAllele2, childAllele1, childAllele2, false);
};

//! Adds a case father duo.
void AllCountedData::addCaseFatherDuo(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * father, Subject * child)
{
  caseFatherDuos.addDuo(snpId, reverseAlleleLabels, father, child);
};

//! Adds a case/mother duo for phased data (using SHAPEIT).
void AllCountedData::addCaseMotherDuoPhased(const unsigned int & snpId, const bool & reverseAlleleLabels, bool & motherAllele1, bool & motherAllele2, bool & childAllele1, bool & childAllele2)
{
	caseMotherDuos.addDuoPhased(snpId, reverseAlleleLabels, motherAllele1, motherAllele2, childAllele1, childAllele2, true);
};

//! Adds a case mother duo.
void AllCountedData::addCaseMotherDuo(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * mother, Subject * child)
{
  caseMotherDuos.addDuo(snpId, reverseAlleleLabels, mother, child);
};

//! Adds a control father duo.
void AllCountedData::addControlFatherDuo(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * father, Subject * child)
{
  controlFatherDuos.addDuo(snpId, reverseAlleleLabels, father, child);
};

//! Adds a control mother duo.
void AllCountedData::addControlMotherDuo(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * mother, Subject * child)
{
  controlMotherDuos.addDuo(snpId, reverseAlleleLabels, mother, child);
};

//! Adds parents of a case.
void AllCountedData::addParentsOfCase(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * father, Subject * mother)
{
  parentsOfCases.addParents(snpId, reverseAlleleLabels, father, mother);
};

//! Adds parents of control.
 void AllCountedData::addParentsOfControl(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * father, Subject * mother)
{
  parentsOfControls.addParents(snpId, reverseAlleleLabels, father, mother);
};

//! Adds father of case.
void AllCountedData::addFatherOfCase(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * father)
{
	fathersOfCases.addSubject(snpId, reverseAlleleLabels, father);
};

//! Adds mother a case.
void AllCountedData::addMotherOfCase(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * mother)
{
  mothersOfCases.addSubject(snpId, reverseAlleleLabels, mother);
};

//! Adds a case.
void AllCountedData::addCase(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * subject)
{
  cases.addSubject(snpId, reverseAlleleLabels, subject);
};

//! Adds a control.
void AllCountedData::addControl(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * subject)
{
  controls.addSubject(snpId, reverseAlleleLabels, subject);
};

//! Outputs to file all counted genotype data to the appropriate files.
void AllCountedData::outputResults(const unsigned int & noOfSnps, const unsigned int & totalNoOfPedigrees, const unsigned int & splitSNPOutput)
{

	caseParentTrios.outputResultsFile(noOfSnps, totalNoOfPedigrees, splitSNPOutput);
	caseMotherDuos.outputResultsFile(noOfSnps, totalNoOfPedigrees, splitSNPOutput);
	caseFatherDuos.outputResultsFile(noOfSnps, totalNoOfPedigrees, splitSNPOutput);
	cases.outputResultsFile(noOfSnps, totalNoOfPedigrees, splitSNPOutput);
	parentsOfCases.outputResultsFile(noOfSnps, totalNoOfPedigrees, splitSNPOutput);
	mothersOfCases.outputResultsFile(noOfSnps, totalNoOfPedigrees, splitSNPOutput);
	fathersOfCases.outputResultsFile(noOfSnps, totalNoOfPedigrees, splitSNPOutput);
	parentsOfControls.outputResultsFile(noOfSnps, totalNoOfPedigrees, splitSNPOutput);
	controlMotherDuos.outputResultsFile(noOfSnps, totalNoOfPedigrees, splitSNPOutput);
	controlFatherDuos.outputResultsFile(noOfSnps, totalNoOfPedigrees, splitSNPOutput);
	controls.outputResultsFile(noOfSnps, totalNoOfPedigrees, splitSNPOutput);
	
};

//! Outputs to screen a summary of the counted genotypes contained in each file.
void AllCountedData::outputFileSummaries(const unsigned int & noOfSnps, const unsigned int & totalNoOfPedigrees)
{

	caseParentTrios.outputFileSummary(noOfSnps, totalNoOfPedigrees);
	caseMotherDuos.outputFileSummary(noOfSnps, totalNoOfPedigrees);
	caseFatherDuos.outputFileSummary(noOfSnps, totalNoOfPedigrees);
	cases.outputFileSummary(noOfSnps, totalNoOfPedigrees);
	parentsOfCases.outputFileSummary(noOfSnps, totalNoOfPedigrees);
	mothersOfCases.outputFileSummary(noOfSnps, totalNoOfPedigrees);
	fathersOfCases.outputFileSummary(noOfSnps, totalNoOfPedigrees);
	parentsOfControls.outputFileSummary(noOfSnps, totalNoOfPedigrees);
	controlMotherDuos.outputFileSummary(noOfSnps, totalNoOfPedigrees);
	controlFatherDuos.outputFileSummary(noOfSnps, totalNoOfPedigrees);
	controls.outputFileSummary(noOfSnps, totalNoOfPedigrees);

};

//! Outputs one line of SNP data to file for the different counted genotype groups.
void AllCountedData::outputResultOneSnp(const unsigned int & snpId, const unsigned int & splitSNPOutput, bool & pooHaps)
{
	
	caseParentTrios.outputResultFileOneSnp(snpId, splitSNPOutput, pooHaps);
	caseMotherDuos.outputResultFileOneSnp(snpId, splitSNPOutput, pooHaps);
	caseFatherDuos.outputResultFileOneSnp(snpId, splitSNPOutput, pooHaps);
	cases.outputResultFileOneSnp(snpId, splitSNPOutput);
	parentsOfCases.outputResultFileOneSnp(snpId, splitSNPOutput);
	mothersOfCases.outputResultFileOneSnp(snpId, splitSNPOutput);
	fathersOfCases.outputResultFileOneSnp(snpId, splitSNPOutput);
	parentsOfControls.outputResultFileOneSnp(snpId, splitSNPOutput);
	controlMotherDuos.outputResultFileOneSnp(snpId, splitSNPOutput);
	controlFatherDuos.outputResultFileOneSnp(snpId, splitSNPOutput);
	controls.outputResultFileOneSnp(snpId, splitSNPOutput);

};

//! Outputs the EMIM parameter file, the first bit.
void AllCountedData::outputParameterFile(ofstream & fileOutPut, const unsigned int & noOfSnps, const bool & useCPG)
{
	
	fileOutPut << "-----------INPUT DATAFILES------------------------------------------------\n"
		<< caseParentTrios.hasData() << "   << caseparenttrios.dat file (0=no, 1=yes, 2=yes, using haplotype estimates)\n"
		<< parentsOfCases.hasData() << "   << caseparents.dat file (0=no, 1=yes)\n"
		<< caseMotherDuos.hasData() << "   << casemotherduos.dat file (0=no, 1=yes, 2=yes, using haplotype estimates)\n"
		<< caseFatherDuos.hasData() << "   << casefatherduos.dat file (0=no, 1=yes, 2=yes, using haplotype estimates)\n"
		<< mothersOfCases.hasData() << "   << casemothers.dat file (0=no, 1=yes)\n"
		<< fathersOfCases.hasData() << "   << casefathers.dat file (0=no, 1=yes)\n"
		<< cases.hasData() << "   << cases.dat file (0=no, 1=yes)\n"
		<< parentsOfControls.hasData() << "   << conparents.dat file (0=no, 1=yes)\n"
		<< controlMotherDuos.hasData() << "   << conmotherduos.dat file (0=no, 1=yes)\n"
		<< controlFatherDuos.hasData() << "   << confatherduos.dat file (0=no, 1=yes)\n"
		<< controls.hasData() << "   << cons.dat file (0=no, 1=yes)\n"
		<< noOfSnps << "   << no of SNPs in each file\n"
		       << "------------------PARAMETER RESTRICTIONS----------------------------------\n"
		<< 0 << "   << fix allele freq A (0=no, 1=yes)\n"
		<< !useCPG << "   << assume HWE and random mating (0=no=estimate 6 mu parameters, 1=yes)\n"
		<< 0 << "   << assume parental allelic exchangeability (0=no, 1=yes)\n"
		<< useCPG << "   << use CPG likelihood (9 mu parameters)\n";

		
};

//! Outputs the EMIM parameter file, the second bit with definable options.
void ProcessData::outputParameterFile(AllCountedData & allCountedData, const unsigned int & noOfSnps) const
{

	string parameterFilename = "emimparams.dat";
	ofstream fileOutput(parameterFilename.c_str());
	
	allCountedData.outputParameterFile(fileOutput, noOfSnps, (probandFileName != ""));

	fileOutput << childGenotype << "   << estimate R1 (0=no, 1=yes)\n"
			<< childGenotype << "   << estimate R2 (0=no, 1=yes)\n"
			<< 0 << "   << R2=R1 (0=no, 1=yes)\n"
			<< childTrend << "   << R2=R1squared	(0=no, 1=yes)\n"
			<< motherGenotype << "   << estimate S1 (0=no, 1=yes)\n"
			<< motherGenotype << "   << estimate S2 (0=no, 1=yes)\n"
			<< 0 << "   << S2=S1 (0=no, 1=yes)\n"
			<< motherTrend << "   << S2=S1squared	(0=no, 1=yes)\n"
			<< imprintingMaternal << "   << estimate Im (0=no, 1=yes)\n"
			<< imprintingPaternal << "   << estimate Ip (0=no, 1=yes)\n"
			<< 0 << "   << estimate gamma11 (0=no, 1=yes)\n"
			<< 0 << "   << estimate gamma12 (0=no, 1=yes)\n"
			<< 0 << "   << estimate gamma21 (0=no, 1=yes)\n"
			<< 0 << "   << estimate gamma22 (0=no, 1=yes)\n"
			<< 0 << "   << gamma22=gamma12=gamma21=gamma11 (0=no, 1=yes)\n"
			<< "---------------OTHER PARAMETERIZATIONS------------------------------------\n"
			<< imprintingMaternalWeinberg << "   << estimate Weinberg (1999b) Im (0=no, 1=yes)\n"
			<< imprintingPaternalWeinberg << "   << estimate Weinberg (1999b) Ip (=Li 2009 Jm) (0=no, 1=yes)\n"
			<< 0 << "   << estimate Sinsheimer (2003) gamma01 (0=no, 1=yes)\n"
			<< 0 << "   << estimate Sinsheimer (2003) gamma21 (0=no, 1=yes)\n"
			<< 0 << "   << estimate Palmer (2006) match parameter (0=no, 1=yes)\n"
			<< 0 << "   << estimate Li (2009) conflict parameter Jc (0=no, 1=yes)\n";

	fileOutput.close();
};

//! Outputs the allele names that were chosen as risk alleles
void ProcessData::outputRiskAlleles(CountAlleles *& countAlleles) const
{
	bool allele1IsRiskAllele;
	map<unsigned int, string>::const_iterator sn;
	string snpName, riskAlleleName;
	map<string, bool> alleleNames; //allele name, and then coding used in PREMIM, e.g SNP 100 may be A -> 0 -> allele "1" and G -> 1 -> allele "2"
	                               // for a binary file these have been set such that, e.g. SNP 100 may be 0 -> allele "2" and 1 -> allele "1"
	map<string, bool>::const_iterator a;
	
	ofstream riskAlleleOutput(riskAlleleOutFileName.c_str());

	//loop thro' all of the SNPs
	for(unsigned int snp = 1; snp <= noOfSnps; ++snp)
	{	
		//Need to check if the SNP has a risk allele specified in risk allele file
		sn = riskAlleleNames.find(snp);

		if(sn != riskAlleleNames.end())
		{
			riskAlleleName = sn->second;
		}
		else
		{
			//Check if the first allele is the risk allele. First allele according to PREMIM, i.e. the first allele encounted in file this may or may not be allele "1" as given in EMIM docs
			allele1IsRiskAllele = (countAlleles->isSnpAllele1Minor(snp) && !useMajorAlleleAsRisk) || (!countAlleles->isSnpAllele1Minor(snp) && useMajorAlleleAsRisk);

			//get the two Allele names for the SNP
			alleleNames = snpAlleleIds.getSnpAlleleNames(snp);
		
			a = alleleNames.begin();

			//check if first allele in list is allele "1" and if it is the risk allele
			//or if first allele in list is allele "2" and it is the risk allele
			if((!a->second && allele1IsRiskAllele) || (a->second && !allele1IsRiskAllele)) riskAlleleName = a->first;
			else
			{
				a++;
				if(a != alleleNames.end()) riskAlleleName = a->first;
				else riskAlleleName = "Unknown";
			};
		};

		snpName = snpNames.find(snp)->second; //rs... number name for SNP

		riskAlleleOutput << snp << "\t" << snpName << "\t" << riskAlleleName << "\n";
	};

	riskAlleleOutput.close();
};

//!Set the list of risk alleles to be used
void ProcessData::setRiskAlleles()
{
	ifstream readRiskAlleleFile(riskAlleleInFileName.c_str());

	if(!readRiskAlleleFile.is_open())
	{
		outErr("\nCannot open risk allele file "); outErr(riskAlleleInFileName); outErr("!\n");
		exit(1);
	};

	string snpNoName, snpAlleleName, alleleName;
	string prevSnpNoName = "";
	unsigned int snpNo = 1;

	do{

		readRiskAlleleFile >> snpNoName >> snpAlleleName >> alleleName;

		if(prevSnpNoName == snpNoName) break;

		riskAlleleNames[snpNo] = alleleName;

		prevSnpNoName = snpNoName;
		snpNo++;

	}while(!readRiskAlleleFile.eof());

	if(riskAlleleNames.size() != noOfSnps)
	{
		outErr("Incorrect number of SNPs in the risk allele file!\n");
		exit(1);
	};

	readRiskAlleleFile.close();
};

//!Returns whether the alleles for a SNP must have their labelling reversed "1" <--> "2" 
bool ProcessData::getReverseAlleles(CountAlleles *& countAlleles, const unsigned int & snpNo) const
{
	if(riskAlleleInFileName == "")
	{
		return (countAlleles->isSnpAllele1Minor(snpNo) && !useMajorAlleleAsRisk) || (!countAlleles->isSnpAllele1Minor(snpNo) && useMajorAlleleAsRisk);
	};

	map<unsigned int, string>::const_iterator ra = riskAlleleNames.find(snpNo);
	if(ra == riskAlleleNames.end())
	{
		out("Risk allele not found for SNP number "); out(snpNo); out("!\n");
		exit(1);
	};

	map<string, bool> snpAls = snpAlleleIds.getSnpAlleleNames(snpNo);
	string firstAl = snpAls.begin()->first;
	bool isMinor = snpAls.begin()->second;

	//check if alleles need to be swapped according to which allele is named in the risk allele file and which is already chosen as the risk allele
	if((ra->second == firstAl && isMinor) || (ra->second != firstAl && !isMinor)) return false;

	return true;
};

//! Adds a subject to the pedigree.
void AllPedigreeData::addSubjectToPedigree(unsigned int & subjectId, Subject * subject, unsigned int & pedigreeId)
{
	map<unsigned int, Pedigree *>::iterator p = pedigrees.find(pedigreeId);

	//if pedigree exists add the subject to it otherwise add subject to a new pedigree
	if(p != pedigrees.end())
	{
		p->second->addSubject(subjectId, subject);
	}
	else
	{
		Pedigree * pedigree = new Pedigree(pedigreeId);
		pedigree->addSubject(subjectId, subject);
		pedigrees[pedigreeId] = pedigree;
	};

};

//! Outputs basic stats of the pedigree to screen.
void AllPedigreeData::outputSummary()
{
	double noOfPedigrees = (double)pedigrees.size();
	double meanPedigreeSize = 0;
	double standardDeviationPedigreeSize = 0;

	for(map<unsigned int, Pedigree *>::const_iterator ped = pedigrees.begin(); ped != pedigrees.end(); ++ped)
	{
		meanPedigreeSize += ped->second->getNumberOfSubjects();
	};
	meanPedigreeSize /= noOfPedigrees;

	for(map<unsigned int, Pedigree *>::const_iterator ped = pedigrees.begin(); ped != pedigrees.end(); ++ped)
	{
		standardDeviationPedigreeSize += (ped->second->getNumberOfSubjects()- meanPedigreeSize)*(ped->second->getNumberOfSubjects()- meanPedigreeSize);
	};
	standardDeviationPedigreeSize /= (noOfPedigrees - 1);
	standardDeviationPedigreeSize = sqrt(standardDeviationPedigreeSize);

	out("Number of pedigrees: "); out(noOfPedigrees); out("\n");
	out("Mean pedigree size: "); out(meanPedigreeSize); out("\n");
	out("Standard deviation of pedigree size: "); out(standardDeviationPedigreeSize); out("\n\n");
};

//! Unmarks subjects in pedigree, used when allowing extras trios.
void AllPedigreeData::restoreAllSubjectsToAllPedigrees()
{
	for(map<unsigned int, Pedigree *>::iterator p = pedigrees.begin(); p != pedigrees.end(); ++p)
	{
		p->second->restoreAllSubjectsToPedigree();
	};

};

//! checks if a subject exists and if it does whether it is affected or not
void AllPedigreeData::subjectExistsAndAffected(const unsigned int & pedId, const unsigned int & subjectId, MapIds & subjectIds) 
{
	bool exists = false;
	bool affected = false;
	Subject * subject;
	
	for(map<unsigned int, Pedigree *>::iterator p = pedigrees.begin(); p != pedigrees.end(); ++p)
	{
		subject = p->second->getSubject(subjectId, exists);
		if(exists)
		{
			affected = subject->getAffected();
			break;
		};
	};

	if(!exists || !affected)
	{
		string subjectName = subjectIds.getName(subjectId);
		
		if(!exists) {outErr("Warning: proband subject "); outErr(subjectName); outErr(" does not exist!\n");}
		else if(!affected) {outErr("Warning: proband subject "); outErr(subjectName); outErr(" is not affected!\n");};
	};
	
};

//! Adds a proband subject to counted data as a trio or something else searching for the best subset in the pedigree
bool ProcessData::addProcessedPedigreeProbandSubject(bool & found, bool & foundTrio, bool & foundParentsCon,
	unsigned int & probandSubjectId, const unsigned int & snpId, const bool & reverseAlleleLabels, Pedigree * pedigree, AllCountedData & allCountedData)
{
	//find the proband subject in the pedigree 
	Subject * probandSubject = pedigree->getSubject(probandSubjectId);
	if(probandSubjectId == 0) return false;

	//try to find a case parent trio
	Trio aTrio;
	pedigree->findCaseParentTrio(aTrio, snpId, found, probandSubject);
	if(found)
	{
		allCountedData.addCaseParentTrio(snpId, reverseAlleleLabels, aTrio.father, aTrio.mother, aTrio.child);
		foundTrio = true;

		if(extraAffectedTrios)
		{			
			aTrio.child->removeFromPedigree();
		};
	};

	if(found) return true;

	//try to find a case mother duo
	pair<Subject *, Subject *> caseMotherDuo;
	pedigree->findCaseMotherDuo(caseMotherDuo, snpId, found, probandSubject);
	if(found)
	{
		allCountedData.addCaseMotherDuo(snpId, reverseAlleleLabels, caseMotherDuo.second, caseMotherDuo.first);
		return true;
	};

	//try to find a case father duo
	pair<Subject *, Subject *> caseFatherDuo;
	pedigree->findCaseFatherDuo(caseFatherDuo, snpId, found, probandSubject);
	if(found)
	{
		allCountedData.addCaseFatherDuo(snpId, reverseAlleleLabels, caseFatherDuo.second, caseFatherDuo.first);
		return true;
	};

	//try to find a case
	Subject * aCase = 0;
	pedigree->findCase(aCase, snpId, found, probandSubject);
	if(found)
	{
		allCountedData.addCase(snpId, reverseAlleleLabels, aCase);
		return true;
	};

	//try to find parents of a case
	pair<Subject *, Subject *> parentsOfCase;
	pedigree->findParentsOfCase(parentsOfCase, snpId, found, probandSubject);
	if(found)
	{
		allCountedData.addParentsOfCase(snpId, reverseAlleleLabels, parentsOfCase.first, parentsOfCase.second);
		return true;
	};

	//try to find a mother of a case
	Subject * motherOfCase = 0;
	pedigree->findMotherOfCase(motherOfCase, snpId, found, probandSubject);
	if(found)
	{
		allCountedData.addMotherOfCase(snpId, reverseAlleleLabels, motherOfCase);
		return true;
	};

	//try to find a father of a case
	Subject * fatherOfCase = 0;
	pedigree->findFatherOfCase(fatherOfCase, snpId, found, probandSubject);
	if(found)
	{
		allCountedData.addFatherOfCase(snpId, reverseAlleleLabels, fatherOfCase);
		return true;
	};

	//the proband subject should be affected so we sould not need the controls below, but they are included just in case
	
	//try to find parents of a control
	pedigree->findParentsOfControl(aTrio, snpId, found, probandSubject);
	if(found)
	{
		allCountedData.addParentsOfControl(snpId, reverseAlleleLabels, aTrio.father, aTrio.mother);
		foundParentsCon = true;

		if(extraUnaffectedTrios)
		{			
			aTrio.father->removeFromPedigree();
			aTrio.mother->removeFromPedigree();
			aTrio.child->removeFromPedigree();
		};
	};
	
	if(found) return true;

	//try to find control mother duo
	pair<Subject *, Subject *> controlMotherDuo;
	pedigree->findControlMotherDuo(controlMotherDuo, snpId, found, probandSubject);
	if(found)
	{
		allCountedData.addControlMotherDuo(snpId, reverseAlleleLabels, controlMotherDuo.second, controlMotherDuo.first);
		return true;
	};

	//try to find control father duo
	pair<Subject *, Subject *> controlFatherDuo;
	pedigree->findControlFatherDuo(controlFatherDuo, snpId, found, probandSubject);
	if(found)
	{
		allCountedData.addControlFatherDuo(snpId, reverseAlleleLabels, controlFatherDuo.second, controlFatherDuo.first);
		return true;
	};

	//try to find control
	Subject * control = 0;
	pedigree->findControl(control, snpId, found, probandSubject);
	if(found)
	{
		allCountedData.addControl(snpId, reverseAlleleLabels, control);
		return true;
	};

	return false;
};

//! Analyses a given pedigree to find the best pedigree subset to count.
bool ProcessData::addProcessedPedigree(const unsigned int & snpId, const bool & reverseAlleleLabels, Pedigree * pedigree, AllCountedData & allCountedData, map<unsigned int, set<unsigned int> > & probandSubjectIds)
{
	bool found = false;
	bool foundTrio = false;
	bool foundParentCons = false;

	//add proband subject first
	map<unsigned int, set<unsigned int> >::const_iterator ps = probandSubjectIds.find(pedigree->getPedId());
	if(ps != probandSubjectIds.end())
	{
		unsigned int probandSubjectId = *(ps->second.begin());
		addProcessedPedigreeProbandSubject(found, foundTrio, foundParentCons, probandSubjectId, snpId, reverseAlleleLabels, pedigree, allCountedData);
	};
	
	//if proband subject has found something then we are done - unless we have chosen to find extras and there may be some to find
	if(found && (!extraAffectedTrios || !foundTrio) &&
		        (!extraUnaffectedTrios || !foundParentCons)) return true;

	bool foundOne;
	Trio aTrio;

	//this check is in case a proband subject yields control parents and we are allowing extra control parents
	if(!(found && extraUnaffectedTrios && foundParentCons))
	{

		//try to find an affected child and parents, if extra trios are allowed, keep searching while removing found trios from search					
		do{
			foundOne = false;
			aTrio = pedigree->findCaseParentTrio(snpId, foundOne);
			if(foundOne)
			{
				allCountedData.addCaseParentTrio(snpId, reverseAlleleLabels, aTrio.father, aTrio.mother, aTrio.child);
				found = true;

				if(extraAffectedTrios)
				{					
					aTrio.child->removeFromPedigree();
				};
			};

		}while(extraAffectedTrios && foundOne);

		if(found) return true;

		//try to find a case mother duo
		pair<Subject *, Subject *> caseMotherDuo =  pedigree->findCaseMotherDuo(snpId, found);
		if(found)
		{
			allCountedData.addCaseMotherDuo(snpId, reverseAlleleLabels, caseMotherDuo.second, caseMotherDuo.first);
			return true;
		};

		//try to find a case father duo
		pair<Subject *, Subject *> caseFatherDuo =  pedigree->findCaseFatherDuo(snpId, found);
		if(found)
		{
			allCountedData.addCaseFatherDuo(snpId, reverseAlleleLabels, caseFatherDuo.second, caseFatherDuo.first);
			return true;
		};

		//try to find a case
		Subject * aCase = pedigree->findCase(snpId, found);
		if(found)
		{
			allCountedData.addCase(snpId, reverseAlleleLabels, aCase);
			return true;
		};

		//try to find parents of a case
		pair<Subject *, Subject *> parentsOfCase = pedigree->findParentsOfCase(snpId, found);
		if(found)
		{
			allCountedData.addParentsOfCase(snpId, reverseAlleleLabels, parentsOfCase.first, parentsOfCase.second);
			return true;
		};

		//try to find mother of a case
		Subject * motherOfCase = pedigree->findMotherOfCase(snpId, found);
		if(found)
		{
			allCountedData.addMotherOfCase(snpId, reverseAlleleLabels, motherOfCase);
			return true;
		};

		//try to find father of a case
		Subject * fatherOfCase = pedigree->findFatherOfCase(snpId, found);
		if(found)
		{
			allCountedData.addFatherOfCase(snpId, reverseAlleleLabels, fatherOfCase);
			return true;
		};

	};
		
	//try to find parents of controls and keep find extra ones if we are allowing for extras
	do{		
		foundOne = false;		
		aTrio = pedigree->findParentsOfControl(snpId, foundOne);
		if(foundOne)
		{			
			allCountedData.addParentsOfControl(snpId, reverseAlleleLabels, aTrio.father, aTrio.mother);
			found = true;

			if(extraUnaffectedTrios)
			{				
				aTrio.father->removeFromPedigree();
				aTrio.mother->removeFromPedigree();
				aTrio.child->removeFromPedigree();				
			};
		};

	}while(extraUnaffectedTrios && foundOne);

	if(found) return true;

	//try to find control mother duo
	pair<Subject *, Subject *> controlMotherDuo = pedigree->findControlMotherDuo(snpId, found);
	if(found)
	{
		allCountedData.addControlMotherDuo(snpId, reverseAlleleLabels, controlMotherDuo.second, controlMotherDuo.first);
		return true;
	};

	//try to find control father duo
	pair<Subject *, Subject *> controlFatherDuo = pedigree->findControlFatherDuo(snpId, found);
	if(found)
	{
		allCountedData.addControlFatherDuo(snpId, reverseAlleleLabels, controlFatherDuo.second, controlFatherDuo.first);
		return true;
	};

	//try to find control
	Subject * control = pedigree->findControl(snpId, found);
	if(found)
	{
		allCountedData.addControl(snpId, reverseAlleleLabels, control);
		return true;
	};

	return false;
};

//! Processes a pedigree for parent-of-origin haplotype trios
bool ProcessData::findPooHapsTrio(const unsigned int & snpId, map<unsigned int, map<unsigned int, unsigned int> > & childrenCountForTrio, const bool & reverseAlleleLabels, Pedigree * pedigree, map<unsigned int, set<unsigned int> > & probandSubjectIds)
{
	bool found = false;
	bool foundTrio = false;
	Trio aTrio;

	//add proband subject first
	map<unsigned int, set<unsigned int> >::const_iterator ps = probandSubjectIds.find(pedigree->getPedId());
	if(ps != probandSubjectIds.end())
	{
		unsigned int probandSubjectId = *(ps->second.begin());
		Subject * probandSubject = pedigree->getSubject(probandSubjectId);
		if(probandSubjectId == 0) return false;

		//try to find a case parent trio		
		pedigree->findCaseParentTrio(aTrio, snpId, found, probandSubject);
		if(found)
		{			
			addToChildrenCountForTrioOrDuo(childrenCountForTrio, pedigree->getPedId(), probandSubjectId);
			return true;
		};

		return false;//do not choose another trio - must stick with the proband
	};
	
	//do every subject, check if child of trio, then pick the one with the least missing data later
	return pedigree->findPooHapCaseParentTrios(snpId, childrenCountForTrio);
};

//! Processes a pedigree for parent-of-origin haplotype mother duos
bool ProcessData::findPooHapsDuo(const unsigned int & snpId, map<unsigned int, map<unsigned int, unsigned int> > & childrenCountForDuo, const bool & reverseAlleleLabels, Pedigree * pedigree, map<unsigned int, set<unsigned int> > & probandSubjectIds, const bool & motherDuo)
{
	bool found = false;
	bool foundDuo= false;
	pair<Subject *, Subject *> aDuo; //child, parent

	//add proband subject first
	map<unsigned int, set<unsigned int> >::const_iterator ps = probandSubjectIds.find(pedigree->getPedId());
	if(ps != probandSubjectIds.end())
	{
		unsigned int probandSubjectId = *(ps->second.begin());
		Subject * probandSubject = pedigree->getSubject(probandSubjectId);
		if(probandSubjectId == 0) return false;

		//try to find a duo	
		if(motherDuo) pedigree->findCaseMotherDuo(aDuo, snpId, found, probandSubject);
		else pedigree->findCaseFatherDuo(aDuo, snpId, found, probandSubject);

		if(found)
		{			
			addToChildrenCountForTrioOrDuo(childrenCountForDuo, pedigree->getPedId(), probandSubjectId);
			return true;
		};

		return false;//do not choose another duo - must stick with the proband when given
	};
	
	if(motherDuo) return pedigree->findPooHapCaseMotherDuos(snpId, childrenCountForDuo);
	else return pedigree->findPooHapCaseFatherDuos(snpId, childrenCountForDuo);
};

//! Returns the sex ID corresponding to the sex ID in file.
unsigned int getSexId(string & sexIdName)
{
	if(sexIdName == "1") return 1;
	else if(sexIdName == "2") return 2;

	return 0;
};

//! Returns the affect status corresponding to the affected value in file.
bool getAffected(string & affectedIdName)
{
	if(affectedIdName == "1") return false;
	else if(affectedIdName == "2") return true;

	return false;
};

//! Outputs the chosen options to add to the EMIM parameter file.
void ProcessData::outputOptionInfo(string & fileName, string & mapFileName)
{	
		out("Input file: "); out(fileName); out("\n");
		if(mapFileName != "") {out("Input map file: "); out(mapFileName); out("\n");};
		if(extraAffectedTrios) out("Allowing extra affected case trios per pedigree.\n");
		if(extraUnaffectedTrios) out("Allowing extra unaffected control matings per pedigree.\n");
		if(childGenotype) out("Child genotype analysis set for parameter file (emimparams.dat).\n");
		if(childTrend) out("Child trend analysis set for parameter file (emimparams.dat).\n");
		if(motherGenotype) out("Mother genotype analysis set for parameter file (emimparams.dat).\n");
		if(motherTrend) out("Mother trend analysis set for parameter file (emimparams.dat).\n");
		if(imprintingMaternal) out("Maternal imprinting analysis set for parameter file (emimparams.dat).\n");
		if(imprintingPaternal) out("Paternal imprinting analysis set for parameter file (emimparams.dat).\n");
		if(imprintingMaternalWeinberg) out("Weinberg maternal imprinting analysis set for parameter file (emimparams.dat).\n");
		if(imprintingPaternalWeinberg) out("Weinberg paternal imprinting analysis set for parameter file (emimparams.dat).\n");
		if(riskAlleleInFileName != "") {out("Using risk alleles from file "); out(riskAlleleInFileName); out(".\n");};
		if(riskAlleleOutFileName != "") {out("Outputting risk alleles to file "); out(riskAlleleOutFileName); out(".\n");};		
		if(useMajorAlleleAsRisk) out("Using the major allele as the risk allele.\n");
		if(estimateAlleleFreq)
		{
			if(splitSNPOutput > 0) out("Estimating allele frequencies and writing to file (emimmarkers*.dat).\n");
			else out("Estimating allele frequencies and writing to file (emimmarkers.dat).\n");
		};
		if(pooHaps)
		{
			out("Estimating parent-of-origin with SHAPEIT2.\n");
			out("SHAPEIT2 command: "); out(shapeitCommand); out("\n");
			out("Maximum missing data permitted for trios and duos: "); out(missThres); out("\n");
		};
};

//! Creates marker files with a given number of SNPs using the given marker file.
void ProcessData::createMarkerFiles(const unsigned int & splitSNPOutput, string & markersFile, string & outputDirectory, string & endName)
{
	out("\nCreating marker files in directory "); out(outputDirectory); out(" using marker file "); out(markersFile); out("...\n");
	ifstream markersIn(markersFile.c_str());

	if(!markersIn.is_open())
	{
		outErr("\nCannot read markers file: "); outErr(markersFile); outErr("!\n");
		exit(1);
	};

	string outputFile = outputDirectory + "emimmarkers"+endName+"1.dat";
	ofstream markersOut(outputFile.c_str());

	string freqName;
	string snpIdName;
	unsigned int snpId = 1;

	//loop thro' markers and write to separate files with splitSNPOutput SNPs each
	do{
		markersIn >> snpIdName >> freqName;
		
		updateMarkersFile(markersOut, snpId, splitSNPOutput, outputDirectory, endName);

		markersOut << snpIdName << "\t" << freqName <<"\n";
		
		++snpId;
	}while(!markersIn.eof());

	markersOut.close();
	markersIn.close();
};

//! Changes the SNP input results file.
bool updateInputFile(ifstream & fileIn, const string & filename, string & outputDirectory, const unsigned int & fileNo)
{
	fileIn.close();
	ostringstream nextFilename;
	nextFilename << outputDirectory << filename<< fileNo << ".out";
	string fname = nextFilename.str();
	fileIn.open(fname.c_str());

	return fileIn.is_open();
};

//! Loops through the different SNP results file and writes to one file.
void writeFile(string & outputDirectory, ifstream & fileIn, ofstream & fileOut, const string & filename)
{
	char returnChar = '\n';
	string aLine;
	unsigned int fileNo = 1;

	//loop thro' all the different SNP results files
	do{
		//only write header line if it is from the first file
		getline(fileIn, aLine, returnChar);
		if(fileNo == 1) fileOut << aLine << "\n";

		//read in every line in file and put into the output file
		do{
			getline(fileIn, aLine, returnChar);

			if(fileIn.eof()) break;

			fileOut << aLine << "\n";
		}while(true);

		//move onto the next file
		++fileNo;		

	}while(updateInputFile(fileIn, filename, outputDirectory, fileNo));

};

//! Concatenates the result files in the given directory into just two files.
void ProcessData::createResultsFiles(string & outputDirectory, string & endName)
{
	out("\nCreating result files emimresults.out and emimsummary.out using files in "); out(outputDirectory); out("...\n");

	string outputResultsFile = "emimresults.out";
	ofstream resultsOut(outputResultsFile.c_str());
	string outputSummaryFile = "emimsummary.out";
	ofstream summaryOut(outputSummaryFile.c_str());

	ifstream resultsIn;
	ifstream summaryIn;
	
	//set up the initial results file
	if(!updateInputFile(resultsIn, "emimresults", outputDirectory, 1))
	{
		outErr("\nCannot read results file: "); outErr(outputDirectory); outErr("emimresults1.out!\n");
		exit(1);
	};

	//set up the initial summary results file
	if(!updateInputFile(summaryIn, "emimsummary", outputDirectory, 1))
	{
		outErr("\nCannot read results file: emimsummary1.out!\n");
		exit(1);
	};

	//create the results file from all the SNP results file
	writeFile(outputDirectory, resultsIn, resultsOut, "emimresults");
	resultsIn.close();
	resultsOut.close();

	//create the summary results file from all the SNP summary results file
	writeFile(outputDirectory, summaryIn, summaryOut, "emimsummary");
	summaryIn.close();
	summaryOut.close();
};

//! Processes the pedigree file - determines pedigree file type for analysis.
void ProcessData::process(string & fileName, string & mapFileName)
{
	//output the chosen options
	outputOptionInfo(fileName, mapFileName);

	unsigned int length = (unsigned)fileName.length();
	string fileExtension = ".ped";
	if(length >= 4) fileExtension = fileName.substr(length-4,4);

	//determine the type of pedigree file
	if(fileExtension[0] == '.' &&
		(fileExtension[1] == 'b' || fileExtension[1] == 'B') &&
		(fileExtension[2] == 'e' || fileExtension[2] == 'E') &&
		(fileExtension[3] == 'd' || fileExtension[3] == 'D'))
	{
		
		if(mapFileName != "")
		{
			outErr("\nA map file may not be used with a binary pedigree file!\n");
			exit(1);
		};

		//open file in binary mode		
		readBinaryGenotypeData.open(fileName.c_str(), ios::binary);
		if(!readBinaryGenotypeData.is_open())
		{
			outErr("\nCannot read binary genotype file: "); outErr(fileName); outErr("!\n");
			exit(1);
		};

		fileType = 2;

		char buffer[3];
		readBinaryGenotypeData.read(buffer, 3);

		//check the plink magic numbers for the file type
		//3rd number indicates format of genotype data, 1 => subjects x SNPs, 0 => SNPs x subjects
		unsigned int magicNumber1 = buffer[0];
		unsigned int magicNumber2 = buffer[1];

		if(magicNumber1 != 108 || magicNumber2 != 27)
		{
			out("Detected an old version .bed file, reading data in SNP-major mode...\n");
			snpMajor = true;
			readBinaryGenotypeData.close();
			readBinaryGenotypeData.open(fileName.c_str(), ios::binary);
		};

		//determine binary file type
		unsigned int mode = buffer[2];
		if(mode == 1)
		{
			snpMajor = true;
		}
		else 
		{
			snpMajor = false;			
		};
	
		bitCount = 9;

		//open corresponding pedigree file for the .bed file
		string famFileName = fileName.substr(0, length-4) + ".fam"; 
		readPedigree.open(famFileName.c_str());

		if(!readPedigree.is_open())
		{
			outErr("\nCannot read corresponding family file: "); outErr(famFileName); outErr("!\n");
			exit(1);
		};
	}
	else if(((fileExtension[0] == 'g' || fileExtension[0] == 'G') &&
		(fileExtension[1] == 'z' || fileExtension[1] == 'Z') &&
		(fileExtension[2] == 'i' || fileExtension[2] == 'I') &&
		(fileExtension[3] == 'p' || fileExtension[3] == 'P'))
		|| 
		(fileExtension[1] == '.' &&
		(fileExtension[2] == 'g' || fileExtension[1] == 'G') &&
		(fileExtension[3] == 'z' || fileExtension[2] == 'Z')))
	{
#ifdef USING_GZIP
		readGzipPedigree.open(fileName.c_str());
		fileType = 3;
#endif

#ifndef USING_GZIP
		outErr("gzipped files are not handled in this version, try a .ped or .bed file instead!\n");
		exit(1);
#endif

	}
	else
	{
		readPedigree.open(fileName.c_str());
		fileType = 1;
	};

	if(!readPedigree.is_open()
#ifdef USING_GZIP
		&& !readGzipPedigree.good()
#endif	
		)
	{
		outErr("\nCannot read pedigree file: "); outErr(fileName); outErr("!\n");
		exit(1);
	};

	//count the number of SNPs and check the chromosome is in 1-22
	setNoOfSnps(fileName, mapFileName);

	//set the risk alleles if set
	if(riskAlleleInFileName != "") setRiskAlleles();

	//run the analysis of the data and output results
	analyseData(fileName);
};

bool isChromosome1to22(string & chromosome)
{
	if(chromosome == "1" || chromosome == "2" || chromosome == "3" 
		|| chromosome == "4" || chromosome == "5" || chromosome == "6" 
		|| chromosome == "7" || chromosome == "8" || chromosome == "9" 
		|| chromosome == "10" || chromosome == "11" || chromosome == "12" 
		|| chromosome == "13" || chromosome == "14" || chromosome == "15" 
		|| chromosome == "16" || chromosome == "17" || chromosome == "18" 
		|| chromosome == "19" || chromosome == "20" || chromosome == "21" || chromosome == "22")
		return true;

	return false;
};

//! Determine how many SNPs there are from the .map or .bim file and check chromosome.
void ProcessData::setNoOfSnps(string & fileName, string & mapFileName)
{
	bool bim = false;
	unsigned int snpCount = 1;
	string chromosome, snpIdentifier, geneticDistance, basePairPosition;
	string alleleName1, alleleName2;
	string prevSnpIdentifier = "";

	if(mapFileName == "")
	{
		unsigned int length = (unsigned)fileName.length();
		if(fileType == 1) mapFileName = fileName.substr(0, length-4) + ".map";
		else if(fileType == 2)
		{
			mapFileName = fileName.substr(0, length-4) + ".bim";
			bim = true;
			bimFileName = mapFileName;
		}
		else
		{
			string letter = fileName.substr(length-1,1);
			if(letter == "p" || letter == "P") mapFileName = fileName.substr(0, length-5) + ".map";
			else mapFileName = fileName.substr(0, length-3) + ".map";
		};
	};

	ifstream readSnpFile;
	readSnpFile.open(mapFileName.c_str());
	if(!readSnpFile.is_open())
	{
		outErr("\nCannot read SNP file: "); outErr(mapFileName); outErr("!\n");
		exit(1);
	};

	do{

		if(bim)
		{
			readSnpFile >> chromosome >> snpIdentifier >> geneticDistance >> basePairPosition
			>>alleleName1 >> alleleName2;

			//record the allele names so that the name of risk alleles can be output if required
			if(riskAlleleOutFileName != "" || riskAlleleInFileName != "")
			{				
				snpAlleleIds.getSnpAlleleId(snpCount, alleleName2); //add second allele, 1 in binary file, as allele "1"
				snpAlleleIds.getSnpAlleleId(snpCount, alleleName1); //add first allele, 0 in binary file, as allele "2"
			};
		}
		else
		{
			readSnpFile >> chromosome >> snpIdentifier >> geneticDistance >> basePairPosition;			
		};
	
		if(snpIdentifier != prevSnpIdentifier) snpNames[snpCount] = snpIdentifier;
		prevSnpIdentifier = snpIdentifier;
		snpCount++;

		if(!isChromosome1to22(chromosome))
		{
			outErr("EMIM is not designed for use with chromosome "); outErr(chromosome); outErr("!\n\n");
			outErr("Please remove all chromosome "); outErr(chromosome); outErr(" SNPs from your input files!\n\n");
			exit(1);
		};
	}while(!readSnpFile.eof());

	noOfSnps = (unsigned)snpNames.size();

	readSnpFile.close();
};

//! Gets a line of information from the pedigree file, .ped or .bim.
void ProcessData::getSubjectNamesFromFile()
{
	if(fileType == 1 || fileType == 2)//.ped or .bim file
	{
		readPedigree >> pedigreeName >> subjectIdName >> fatherIdName >> motherIdName >> sexIdName >> affectedIdName;
	}
#ifdef USING_GZIP
	else//.gzip file
	{
		readGzipPedigree >> pedigreeName >> subjectIdName >> fatherIdName >> motherIdName >> sexIdName >> affectedIdName;
	};
#endif
};

//! Get genotype from file.
void ProcessData::getAlleleNamesFromFile()
{	
#ifdef USING_GZIP
	if(fileType == 3)
	{
		readGzipPedigree >> allele1Name >> allele2Name; 
	}
	else
#endif
		readPedigree >> allele1Name >> allele2Name; 
	
};

//! Adds genotype data to subject.
void ProcessData::addGenotypeDataToSubject(Subject * subject, const unsigned int & snpId, CountAlleles * countAlleles)
{
	if(fileType == 2)
	{
		addGenotypeDataToSubjectBinary(subject, snpId, countAlleles);
	}
	else
	{
		addGenotypeDataToSubjectPed(subject, snpId, countAlleles);
	};

};

//! Adds genotype data to subject for a SNP for non-binary pedigree files.
void ProcessData::addGenotypeDataToSubjectPed(Subject * subject, const unsigned int & snpId, CountAlleles * countAlleles)
{

	//read in alleles from file
	getAlleleNamesFromFile();

	if(readPedigree.eof()
#ifdef USING_GZIP
		|| readGzipPedigree.eof()
#endif		
		)
	{
		outErr("Error: last row of SNPs incomplete!\n");
		exit(1);
	};

	//if genotype data is not defined then do not record it
	if(allele1Name != "0" && allele2Name != "0")
	{
		//look up allele codes
		allele1Id = snpAlleleIds.getSnpAlleleId(snpId, allele1Name);
		allele2Id = snpAlleleIds.getSnpAlleleId(snpId, allele2Name);

		//add genotype data to the subject
		subject->addGenotype(snpId, allele1Id, allele2Id);
		
		//count allele1 to see if it is the minor allele
		unsigned int allele1Count = 0;
		if(!allele1Id) allele1Count++;
		if(!allele2Id) allele1Count++;
		countAlleles->addToAllele1Count(snpId, allele1Count, subject->getAffected());

	};

};

//! Sets the genotype data for the given subject and SNP for a binary pedigree file.
void ProcessData::addGenotypeDataToSubjectBinary(Subject * subject, const unsigned int & snpId, CountAlleles * countAlleles)
{
	
	//read in the next piece of data
	if(bitCount == 9)
	{
		
		readBinaryGenotypeData.read(buffer, 1);
		if(readBinaryGenotypeData.eof())
		{			
			outErr("Error: binary SNP file (.bed) is incomplete!\n");
			exit(1);
		};

		aBit = buffer[0];

		bitCount = 1;
	};

	allele1 = aBit & one; //read the least significant bit				
	aBit = aBit >> 1; //shift bits to the right
	allele2 = aBit & one; //read the new least significant bit				
	aBit = aBit >> 1; //shift bits to the right for next time

	bitCount += 2;	

	if(allele1 == 1 && allele2 == 1)
	{
		allele1 = 0; allele2 = 0;
	}
	else if(allele1 == 0 && allele2 == 0)
	{
		allele1 = 1; allele2 = 1;
	};
		
	//add genotype data to the subject, (switch 0 and 1 so that 1 is most likely to be minor allele)
	subject->addGenotype(snpId, (allele1 == 1), (allele2 == 1));	
	
	if(allele1 != 1 || allele2 != 0)
	{
		//count allele1 to see if it is the minor allele
		unsigned int allele1Count = 2 - allele1 - allele2; 
		countAlleles->addToAllele1Count(snpId, allele1Count, subject->getAffected());
	};
	
};


//add the IDs of the subjects in the proband file first so that they appear
//in the list of subjects first and are then given priority chosen for analysis
void addProbandSubjectIds(string & probandFileName, MapIds & subjectIds, map<unsigned int, set<unsigned int> > & probandSubjectIds)
{
	ifstream probandFile;
	string pedigreeName, subjectIdName, subjectIdName2;
	unsigned int subjectId;
	set<unsigned int> pedSubjectIds;

	probandFile.open(probandFileName.c_str());	
	
	if(!probandFile.is_open())
	{
		outErr("Proband file "); outErr(probandFileName); outErr(" not found!\n");
		exit(1);
	};		

	do{
		probandFile >> pedigreeName >> subjectIdName;

		if(probandFile.eof()) break;

		subjectIdName2 = pedigreeName + "-" + subjectIdName;
		subjectId = subjectIds.getId(subjectIdName2);

		//add subject Id to the coressponding pedigree Id
		map<unsigned int, set<unsigned int> >::iterator psi = probandSubjectIds.find(subjectId);
		if(psi != probandSubjectIds.end())
		{
			psi->second.insert(subjectId);
		}
		else
		{
			pedSubjectIds.clear();
			pedSubjectIds.insert(subjectId);
			probandSubjectIds[subjectId] = pedSubjectIds;//pedigreeName + " " + subjectIdName; 
		};
	}while(!probandFile.eof());
	
	probandFile.close();
};

//! This function contains the main loop for processing the pedigree file
void ProcessData::analyseData(string & fileName)
{
	MapIds subjectIds;
	MapIds pedigreeIds;	
	CountAlleles * countAlleles = new CountAlleles();
	
	pedigreesNotCounted = 0;
	map<unsigned int, unsigned int> noOfpedigreeSubjects; //pedigree id, subjects count

	unsigned int subjectId, fatherId, motherId, pedigreeId, sexId;
	bool affected;
	bool reverseAlleles;

	//open file to record the estimated allele frequencies
	if(estimateAlleleFreq)
	{
		//determine if the output is to be split across different files
		string alFrFilename;
		if(splitSNPOutput == 0) alFrFilename = outputDirectory+"emimmarkers"+endName+".dat";
		else alFrFilename = outputDirectory+"emimmarkers"+endName+"1.dat";

		markersFile.open(alFrFilename.c_str());		
	};	

	if(probandFileName != "") addProbandSubjectIds(probandFileName, subjectIds, probandSubjectIds);

	Subject * subject;
	
	unsigned int snpId = 1;
	unsigned int count = 0;
	unsigned int totalNoOfSubjects = 0;
	unsigned int totalNoOfMales = 0;
	unsigned int totalNoOfFemales = 0;
	unsigned int totalNoOfUnknownSex = 0;
	unsigned int totalNoOfSubjectsAffected = 0;
	
	//read in pedigree data
	do{

		//read in all of the details of a subject and map Ids of parents and the subject
		getSubjectNamesFromFile();

		if(readPedigree.eof()
#ifdef USING_GZIP
			|| readGzipPedigree.eof()
#endif			
			) break; //check if the last row has been past
		
		pedigreeId = pedigreeIds.getId(pedigreeName);

		subjectIdName = pedigreeName + "-" + subjectIdName;
		subjectId = subjectIds.getId(subjectIdName);

		if(fatherIdName == "0") fatherId = 0;
		else
		{
			fatherIdName = pedigreeName + "-" + fatherIdName;
			fatherId = subjectIds.getId(fatherIdName);
		};

		if(motherIdName == "0") motherId = 0;
		else
		{
			motherIdName = pedigreeName + "-" + motherIdName;
			motherId = subjectIds.getId(motherIdName);
		};

		sexId = getSexId(sexIdName);

		affected = getAffected(affectedIdName);
	
		if(affected) totalNoOfSubjectsAffected++;
		
		//create a subject
		if(fileType == 2 && snpMajor)
		{			
			subject = new SubjectOneSnp(subjectId, fatherId, motherId, sexId, affected);			

			orderListOfSubjects.push_back(subject);		
		}		
		else
		{
			subject = new SubjectMultiSnp(subjectId, fatherId, motherId, sexId, affected);
		
			if(fileType == 2) bitCount = 9; //start reading data from the next byte

			//add all SNP data to the subject
			snpId = 1;
			do{
				addGenotypeDataToSubject(subject, snpId, countAlleles);

				snpId++;
			}while(snpId <= noOfSnps);

		};

		//add the subject to the pedigree to which they belong
		allPedigreeData.addSubjectToPedigree(subjectId, subject, pedigreeId);

		totalNoOfSubjects++;
		if(sexId == 1) totalNoOfMales++;
		else if(sexId == 2) totalNoOfFemales++;
		else totalNoOfUnknownSex++;

	}while(!readPedigree.eof()
#ifdef USING_GZIP
		&& !readGzipPedigree.eof()
#endif		
		);

	//check all the proband subjects
	for(map<unsigned int, set<unsigned int> >::const_iterator psids = probandSubjectIds.begin(); psids != probandSubjectIds.end(); ++psids)
	{
		for(set<unsigned int>::const_iterator s = psids->second.begin(); s != psids->second.end(); ++s)
		{
			allPedigreeData.subjectExistsAndAffected(psids->first, *s, subjectIds);
		};
	};

	//deal with trios and duos for poo (parent of origin) with haplotype estimation
	if(pooHaps) processPooHaps(fileName, countAlleles); //must be a binary PLINK file

	if(fileType == 2 && snpMajor)
	{		
		for(unsigned int snpm = 1; snpm <= noOfSnps; ++snpm)
		{	
			bitCount = 9; //start reading data from the next byte since in SNP-major mode for the next SNP

			//update genotype data for given SNP for every subject
			for(list<Subject *>::const_iterator s = orderListOfSubjects.begin(); s != orderListOfSubjects.end(); ++s)
			{
				addGenotypeDataToSubject(*s, snpm, countAlleles);
			};

			//check if allele 1 is minor
			countAlleles->updateIsAllele1Minor();
			reverseAlleles = getReverseAlleles(countAlleles, snpm);
		
			//loop through the pedigrees adding trios, duos etc
			for(map<unsigned int, Pedigree *>::iterator ped = allPedigreeData.pedigrees.begin(); ped != allPedigreeData.pedigrees.end(); ++ped)
			{
				if(!pooHaps)
				{
					if(!(addProcessedPedigree(1, reverseAlleles, ped->second, allCountedData, probandSubjectIds)))
					{
						pedigreesNotCounted++;
					};
				}
				else if(!ped->second->isPooHapProcessed())
				{
					if(!(addProcessedPedigree(snpm, reverseAlleles, ped->second, allCountedData, probandSubjectIds)))
					{
						pedigreesNotCounted++;
					};
				};
			};

			if(!pooHaps) allCountedData.outputResultOneSnp(snpm, splitSNPOutput, pooHaps);

			//reset the allele1 count for the next SNP to see if it is minor
			if(estimateAlleleFreq) countAlleles->outputMarkersFileOneLine(markersFile, snpm, splitSNPOutput, outputDirectory, endName, useMajorAlleleAsRisk, riskAlleleNames, snpAlleleIds, (riskAlleleInFileName != ""));
			countAlleles->clearAllele1Counts();

			allPedigreeData.restoreAllSubjectsToAllPedigrees();			
		};
	}
	else //all data was added for .ped file
	{	
		//check if allele 1 is minor for each SNP
		countAlleles->updateIsAllele1Minor();

		//process pedigrees
		for(map<unsigned int, Pedigree *>::iterator ped = allPedigreeData.pedigrees.begin(); ped != allPedigreeData.pedigrees.end(); ++ped)
		{			
			//loop thro' all of the SNPs
			for(unsigned int snp = 1; snp <= noOfSnps; ++snp)
			{				
				reverseAlleles = getReverseAlleles(countAlleles, snp);
				if(!(addProcessedPedigree(snp, reverseAlleles, ped->second, allCountedData, probandSubjectIds))) pedigreesNotCounted++;

				//all subjects may be used for genotype groups for the next SNP
				ped->second->restoreAllSubjectsToPedigree();
			};						
		};	
	};

	//output list of allele names for each SNP that is the used for the risk allele 
	if(riskAlleleOutFileName != "") outputRiskAlleles(countAlleles);

	double femalePercent, malePercent, unknownSexPercent, affectedPercent;
	malePercent = (((double)(totalNoOfMales))/((double)(totalNoOfSubjects)))*100;
	femalePercent = (((double)(totalNoOfFemales))/((double)(totalNoOfSubjects)))*100;
	unknownSexPercent = (((double)(totalNoOfUnknownSex))/((double)(totalNoOfSubjects)))*100;
	affectedPercent = (((double)(totalNoOfSubjectsAffected))/((double)(totalNoOfSubjects)))*100;

	out("\nNumber of subjects: "); out(totalNoOfSubjects);
	out("\n          Males: "); out(totalNoOfMales); out(" ("); out(malePercent); out("%)");
	out("\n          Females: "); out(totalNoOfFemales); out(" ("); out(femalePercent); out("%)");
	out("\n          Unknown sex: "); out(totalNoOfUnknownSex); out(" ("); out(unknownSexPercent); out("%)");
	out("\n          Affected: "); out(totalNoOfSubjectsAffected); out(" ("); out(affectedPercent); out("%)");
	out("\n          Unaffected: "); out(totalNoOfSubjects-totalNoOfSubjectsAffected); out(" ("); out(100-affectedPercent); out("%)");
	out("\nNumber of SNPs: "); out(noOfSnps); out("\n\n");
	
	allPedigreeData.outputSummary();

	if(estimateAlleleFreq)
	{
		if(fileType == 2 && snpMajor) markersFile.close();
		else
		{
			countAlleles->outputMarkersFile(markersFile, noOfSnps, splitSNPOutput, outputDirectory, endName, useMajorAlleleAsRisk, riskAlleleNames, snpAlleleIds, (riskAlleleInFileName != ""));
		};

		out("File name: emimmarkers.dat\n");
		out("Number of allele SNP frequencies estimated: "); out(noOfSnps); out("\n\n");
	};

	unsigned int totalNoOfPedigrees = (unsigned)allPedigreeData.pedigrees.size();

	if(fileType != 2 || !snpMajor || pooHaps) allCountedData.outputResults(noOfSnps, totalNoOfPedigrees, splitSNPOutput);
	else allCountedData.outputFileSummaries(noOfSnps, totalNoOfPedigrees);

	out("Number of uncounted groups: "); out(pedigreesNotCounted); out("\n\n");
	outputParameterFile(allCountedData, noOfSnps);

	if(fileType == 1 || fileType == 2) readPedigree.close();
	if(fileType == 2) readBinaryGenotypeData.close();
#ifdef USING_GZIP
	if(fileType == 3) readGzipPedigree.close();
#endif

	delete countAlleles;
};

//! Write .bim or .map file for use with SHAPEIT, ecludes missing SNPs
void ProcessData::outputPooHapMapOrBimFile(string & fileName, string & inputFilename, set<unsigned int> & includedSNPs)
{
	bool bim = false;
	unsigned int snpCount = 1;
	string chromosome, snpIdentifier, geneticDistance, basePairPosition, mapFileName, outputMapFilename;
	string alleleName1 = "1", alleleName2 = "2";
	string prevSnpIdentifier = "";
	set<unsigned int>::const_iterator is;

	unsigned int length = (unsigned)fileName.length();
	if(fileType == 1)
	{
		mapFileName = fileName.substr(0, length-4) + ".map";		
	}
	else if(fileType == 2)
	{
		mapFileName = fileName.substr(0, length-4) + ".bim";			
	}
	else
	{
		string letter = fileName.substr(length-1,1);
		if(letter == "p" || letter == "P") mapFileName = fileName.substr(0, length-5) + ".map";
		else mapFileName = fileName.substr(0, length-3) + ".map";		
	};
	
	outputMapFilename = inputFilename  + ".bim";

	ifstream readSnpFile;
	readSnpFile.open(mapFileName.c_str());
	if(!readSnpFile.is_open())
	{
		outErr("\nCannot read SNP desc. file: "); outErr(mapFileName); outErr("!\n");
		exit(1);
	};

	ofstream writeSnpFile(outputMapFilename.c_str());
	unsigned int snpNo = 1;

	do{

		if(fileType == 2)
		{
			readSnpFile >> chromosome >> snpIdentifier >> geneticDistance >> basePairPosition
			>>alleleName1 >> alleleName2;

			////record the allele names so that the name of risk alleles can be output if required
			//if(riskAlleleOutFileName != "" || riskAlleleInFileName != "")
			//{				
			//	snpAlleleIds.getSnpAlleleId(snpCount, alleleName2); //add second allele, 1 in binary file, as allele "1"
			//	snpAlleleIds.getSnpAlleleId(snpCount, alleleName1); //add first allele, 0 in binary file, as allele "2"
			//};
		}
		else
		{
			readSnpFile >> chromosome >> snpIdentifier >> geneticDistance >> basePairPosition;			
		};
		
		is = includedSNPs.find(snpNo);
		if(is != includedSNPs.end())
		{		
			writeSnpFile << chromosome << " " << snpIdentifier << " " << geneticDistance << " " << basePairPosition << " " << alleleName1 << " " << alleleName2 << "\n";
		};

		++snpNo;
	}while(!readSnpFile.eof());

	writeSnpFile.close();
	readSnpFile.close();
};


//! Adds trios using phased haplotype info using SHAPEIT
void ProcessData::processPooHaps(string & fileName, CountAlleles * countAlleles)
{
	if(fileType != 2 || !snpMajor)
	{
		outErr("\nA valid PLINK binary pedigree file (.bed/.bim/.fam) must be used for parent-of-origin effects with SHAPEIT2!\n");
		exit(1);
	};

	double percentThresPoo = 1.0 - missThres;//percentage of the same trios (or duos) chosen from a pedigree (from all SNPs) needed to make a poo trio (or duo) 

	bool reverseAlleles;
	string inputFilename = outputDirectory + "tempForSHAPEIT" + endName;

	map<unsigned int, unsigned int> pooHapTrios; //pedigree ID, child subject ID (best found)
	map<unsigned int, unsigned int> pooHapMotherCaseDuos; //pedigree ID, child subject ID (best found)
	map<unsigned int, unsigned int> pooHapFatherCaseDuos; //pedigree ID, child subject ID (best found)
	map<unsigned int, map<unsigned int, unsigned int> > childrenCountForTrio; //pedigree ID, child subject ID, count -- for a pedigree count of each trio given by child subject ID
	map<unsigned int, map<unsigned int, unsigned int> > childrenCountForMotherCaseDuo; //pedigree ID, child subject ID, count -- for a pedigree count of each mother/case duo given by child subject ID
	map<unsigned int, map<unsigned int, unsigned int> > childrenCountForFatherCaseDuo; //pedigree ID, child subject ID, count -- for a pedigree count of each father/case duo given by child subject ID
	set<unsigned int> includedSNPs;
	map<unsigned int, map<unsigned int, unsigned int> >::const_iterator currentCount;

	char buffer[3];
	double percent;
	
	unsigned int snp = 1;

	//loop thro' all of the SNPs
	snp = 1;
	while(snp <= noOfSnps)
	{	
	
		//add SNP data, adding SNP data one at a time	
		bitCount = 9; //start reading data from the next byte since in SNP-major mode for the next SNP

		//update genotype data for given SNP for every subject
		for(list<Subject *>::const_iterator s = orderListOfSubjects.begin(); s != orderListOfSubjects.end(); ++s)
		{
			addGenotypeDataToSubject(*s, snp, countAlleles);
		};

		//check if allele 1 is minor
		countAlleles->updateIsAllele1Minor(snp);			
		reverseAlleles = getReverseAlleles(countAlleles, snp);

		if(!countAlleles->isSNPDataMissing(snp)) includedSNPs.insert(snp);
			
		//process pedigrees
		for(map<unsigned int, Pedigree *>::iterator ped = allPedigreeData.pedigrees.begin(); ped != allPedigreeData.pedigrees.end(); ++ped)
		{
		
			if(ped->second->getNumberOfSubjects() >= 2)
			{
				
				findPooHapsTrio(snp, childrenCountForTrio, reverseAlleles, ped->second, probandSubjectIds);	

				//all subjects may be used for genotype groups for the mother case duos, although may not be used if used for trios
				ped->second->restoreAllSubjectsToPedigree();

				//find mother case duos
				findPooHapsDuo(snp, childrenCountForMotherCaseDuo, reverseAlleles, ped->second, probandSubjectIds, true);

				//all subjects may be used for genotype groups for the father case duos, although may not be used if used for mother/case duos or trios
				ped->second->restoreAllSubjectsToPedigree();

				//find father case duos
				findPooHapsDuo(snp, childrenCountForFatherCaseDuo, reverseAlleles, ped->second, probandSubjectIds, false);

				//all subjects may be used for genotype groups for the next SNP
				ped->second->restoreAllSubjectsToPedigree();
	
			};//end of skipping ped with 1 subject

			pooHapTrios[ped->first] = 0;
			pooHapMotherCaseDuos[ped->first] = 0;
			pooHapFatherCaseDuos[ped->first] = 0;	
		};//end of loop thro' pedigrees
	
	++snp;
	};//end of loop thro' SNPs

	readBinaryGenotypeData.close(); //close SNP .bed file

	outputPooHapMapOrBimFile(fileName, inputFilename, includedSNPs);

	unsigned int bestChildID = 0;
	unsigned int highestNoSNPs = 0;

	//set up which trios to use for each pedigree, 0 = no trio
	//pedigree ID, child subject ID, count 
	for(map<unsigned int, map<unsigned int, unsigned int> >::const_iterator cct = childrenCountForTrio.begin(); cct != childrenCountForTrio.end(); ++cct)
	{
		bestChildID = 0;
		highestNoSNPs = 0;

		//find trio with least missing data
		for(map<unsigned int, unsigned int>::const_iterator c = cct->second.begin(); c != cct->second.end(); ++c)
		{
			//new best trio
			if(c->second > highestNoSNPs)
			{
				bestChildID = c->first;
				highestNoSNPs = c->second;
			};
		};

		if(bestChildID != 0)
		{
			//add if above missing threshold
			percent = ((double)(highestNoSNPs))/((double)(noOfSnps));

			if(percent > percentThresPoo) //enough trios of this sort for a poo trio
			{
				pooHapTrios[cct->first] = bestChildID;				
				allPedigreeData.pedigrees.find(cct->first)->second->setPooHapProcessed(true);//do not process this pedigree in the normal way							
			};
		};

	};


	//set up which mother case duos to use for each pedigree, 0 = no trio or mother case duo 
	for(map<unsigned int, map<unsigned int, unsigned int> >::const_iterator ccmcd = childrenCountForMotherCaseDuo.begin(); ccmcd != childrenCountForMotherCaseDuo.end(); ++ccmcd)
	{
		bestChildID = 0;
		highestNoSNPs = 0;

		//find mother duo with least missing data
		for(map<unsigned int, unsigned int>::const_iterator c = ccmcd->second.begin(); c != ccmcd->second.end(); ++c)
		{
			//new mother duo
			if(c->second > highestNoSNPs)
			{
				bestChildID = c->first;
				highestNoSNPs = c->second;
			};
		};

		if(bestChildID != 0)
		{
			//add if above missing threshold
			percent = ((double)(highestNoSNPs))/((double)(noOfSnps));

			if(percent > percentThresPoo && pooHapTrios[ccmcd->first] == 0) //enough mother duos of this sort for a poo trio  and not already used for a trio
			{
				pooHapMotherCaseDuos[ccmcd->first] = bestChildID;				
				allPedigreeData.pedigrees.find(ccmcd->first)->second->setPooHapProcessed(true);//do not process this pedigree in the normal way							
			};
		};

	};

	//set up which father case duos to use for each pedigree, 0 = no trio or mother case duo 
	for(map<unsigned int, map<unsigned int, unsigned int> >::const_iterator ccfcd = childrenCountForFatherCaseDuo.begin(); ccfcd != childrenCountForFatherCaseDuo.end(); ++ccfcd)
	{
		bestChildID = 0;
		highestNoSNPs = 0;

		//find mother duo with least missing data
		for(map<unsigned int, unsigned int>::const_iterator c = ccfcd->second.begin(); c != ccfcd->second.end(); ++c)
		{
			//new father duo
			if(c->second > highestNoSNPs)
			{
				bestChildID = c->first;
				highestNoSNPs = c->second;
			};
		};

		if(bestChildID != 0)
		{
			//add if above missing threshold
			percent = ((double)(highestNoSNPs))/((double)(noOfSnps));

			if(percent > percentThresPoo && pooHapTrios[ccfcd->first] == 0 && pooHapMotherCaseDuos[ccfcd->first] == 0) //enough father duos of this sort for a poo trio and not already used for a trio or mother/case duo
			{
				pooHapFatherCaseDuos[ccfcd->first] = bestChildID;				
				allPedigreeData.pedigrees.find(ccfcd->first)->second->setPooHapProcessed(true);//do not process this pedigree in the normal way							
			};
		};

	};

	string outputFamFilename = inputFilename + ".fam";
	string outputBedFilename = inputFilename + ".bed";

	ofstream outFamPooFile(outputFamFilename.c_str());
	ofstream outBedPooFile(outputBedFilename.c_str(), ios::binary);

	//write out initial binary pedigree file bytes, first 2 byte are magic numbers the third to indicate SNP major (subjects x SNPs)
	//char buffer[3];	
	buffer[0] = 108;
	buffer[1] = 27;
	buffer[2] = 1;
	outBedPooFile.write(buffer, 3);

	map<unsigned int, unsigned int>::const_iterator pot = pooHapTrios.begin();
	unsigned int noFoundTrios = 0;

	//create data files to run with SHAPEIT to estimate the haplotypes for the poo trios and duos
	//output trio fam file details
	for(map<unsigned int, Pedigree *>::iterator p = allPedigreeData.pedigrees.begin(); p != allPedigreeData.pedigrees.end(); ++p, ++pot)
	{
		if(pot->second != 0)
		{
			outputPooHapTrioFamFile(p->first, p->second, pot->second, outFamPooFile);		
			++noFoundTrios;
		};
	};

	//output mother/case duo fam file details
	map<unsigned int, unsigned int>::const_iterator pomcd = pooHapMotherCaseDuos.begin();
	unsigned int noFoundMotherDuos = 0;
	//create data files to run with SHAPEIT to estimate the haplotypes for the poo mother duos
	for(map<unsigned int, Pedigree *>::iterator pm = allPedigreeData.pedigrees.begin(); pm != allPedigreeData.pedigrees.end(); ++pm, ++pomcd)
	{
		if(pomcd->second != 0)
		{
			outputPooHapMotherDuoFamFile(pm->first, pm->second, pomcd->second, outFamPooFile);		
			++noFoundMotherDuos;
		};
	};

	//output father/case duo fam file details
	map<unsigned int, unsigned int>::const_iterator pofcd = pooHapFatherCaseDuos.begin();
	unsigned int noFoundFatherDuos = 0;
	//create data files to run with SHAPEIT to estimate the haplotypes for the poo father duos
	for(map<unsigned int, Pedigree *>::iterator pf = allPedigreeData.pedigrees.begin(); pf != allPedigreeData.pedigrees.end(); ++pf, ++pofcd)
	{
		if(pofcd->second != 0)
		{
			outputPooHapFatherDuoFamFile(pf->first, pf->second, pofcd->second, outFamPooFile);		
			++noFoundFatherDuos;
		};
	};

	outputPooHapBedFile(fileName, countAlleles, pooHapTrios, pooHapMotherCaseDuos, pooHapFatherCaseDuos, allPedigreeData, includedSNPs, outBedPooFile);

	outBedPooFile.close();
	outFamPooFile.close();

	if(noFoundTrios > 0 || noFoundMotherDuos > 0 || noFoundFatherDuos > 0)
	{

		//run SHAPEIT on created data files
		string outputFilename = outputDirectory + "tempOutSHAPEIT" + endName;
		string inputMap = "";

		string shapeitCommandGraph = shapeitCommand
			+ " --input-bed " + inputFilename + ".bed " + inputFilename + ".bim " + inputFilename + ".fam";
		if(inputMap != "")	shapeitCommandGraph = shapeitCommandGraph	+ " --input-map " + inputMap;

		shapeitCommandGraph = shapeitCommandGraph + " --output-graph " + outputFilename + ".hgraph"
			+ " --output-log " + outputFilename + "-hgraph.log" 
			+ threadOptions + hapMCMCOptions + hapModelOptions;
		if(otherShapeitOtherOptions != "") shapeitCommandGraph = shapeitCommandGraph + " " + otherShapeitOtherOptions;
		shapeitCommandGraph = shapeitCommandGraph + " >/dev/null 2>&1";

		out("\nCalculating haplotype graph using SHAPEIT2 command:\n"); out(shapeitCommandGraph); out("\n\n");

		unsigned int errorCode = system(shapeitCommandGraph.c_str());

		if(errorCode != 0)
		{
			outErr("Problem executing SHAPEIT2 (error code "); outErr(errorCode); outErr(") when calculating the haplotype graph using:\n\n ");
			outErr(shapeitCommandGraph); outErr("\n\n");
			outErr("Log file: "); outErr(outputFilename); outErr("-hgraph.log:\n\n");
			string shapeitLog = "more " + outputFilename + "-hgraph.log";
			system(shapeitLog.c_str());

			outErr("\nCHECK SHAPEIT2 COMMAND IS SETUP CORRECTLY WITH \"-shapeit\" OPTION, type: \""); outErr(shapeitCommand); outErr("\" to check it runs.\n\n");
			outErr("You may need to download and install SHAPEIT2:\n https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#download\n\n");

			exit(errorCode);		
		};

		string shapeitCommandEst = shapeitCommand + " -convert --input-graph " + outputFilename + ".hgraph";

		if(noSamples == 1) shapeitCommandEst = shapeitCommandEst + " --output-max ";
		else shapeitCommandEst = shapeitCommandEst + " --output-sample ";
			
		shapeitCommandEst = shapeitCommandEst + outputFilename + ".haps " + outputFilename + ".sample"
			+ " --output-log " + outputFilename + "-max.log" 
			+ " " + threadOptions
			+ " >/dev/null 2>&1";

		string hapsFilename = outputFilename + ".haps";
		ifstream readHapsFile;
		string chr, snpName, snpPos, allele1, allele2;
		string  snpNamePrev = "";
		bool fatherAllele1, fatherAllele2, motherAllele1, motherAllele2, childAllele1, childAllele2; 

		out("Calculating haplotype estimates using SHAPEIT2 command:\n"); out(shapeitCommandEst); out("\n\n");

		allCountedData.setTotalSamples(noSamples);

		for(unsigned int sampleNo = 1; sampleNo <= noSamples; ++sampleNo)
		{
			errorCode = system(shapeitCommandEst.c_str());

			if(errorCode != 0)
			{
				outErr("\nProblem executing SHAPEIT2 (error code "); outErr(errorCode); outErr(") when calculating the best estimate haplotypes using ");
				outErr(shapeitCommandEst); outErr("!\n\n");
				outErr("Log file: "); outErr(outputFilename); outErr("-max.log:\n\n");
				string shapeitLog = "more " + outputFilename + "-max.log";
				system(shapeitLog.c_str());
				exit(errorCode);		
			};

			//add trios using phase info		
			readHapsFile.open(hapsFilename.c_str());

			if(!readHapsFile.is_open())
			{
				outErr("\nCannot open haps file from SHAPEIT2: "); outErr(hapsFilename); outErr("!\n");
				exit(1);
			};
		
			snpNamePrev = "";
		
			for(list<unsigned int>::const_iterator snpNo = includedSNPsInShapeIT.begin(); snpNo != includedSNPsInShapeIT.end(); ++snpNo)
			{	
				//read in snp info
				readHapsFile >> chr >> snpName >> snpPos >> allele1 >> allele2; 

				if(snpName == snpNamePrev)
				{
					outErr("Problem reading .haps file: "); outErr(hapsFilename); outErr(" at the SNP after: "); outErr(snpName); outErr("!\n");
					exit(1);
				};

				reverseAlleles = getReverseAlleles(countAlleles, *snpNo);
				countAlleles->addReverseAlleles(*snpNo, reverseAlleles); //keep track of whether alleles were reversed or not for when adjusting duos

				//loop thro' trios
				for(unsigned int trioNo = 1; trioNo <= noFoundTrios; ++trioNo)
				{
					readHapsFile >> fatherAllele1 >> fatherAllele2 >> motherAllele1 >> motherAllele2 >> childAllele1 >> childAllele2; 

					allCountedData.addCaseParentTrioPhased(*snpNo, reverseAlleles, fatherAllele1, fatherAllele2, motherAllele1, motherAllele2, childAllele1, childAllele2); //will divide by total later to get ave
				};

				//loop thro' case/mother duos
				for(unsigned int motherDuoNo = 1; motherDuoNo <= noFoundMotherDuos; ++motherDuoNo)
				{
					readHapsFile >> motherAllele1 >> motherAllele2 >> childAllele1 >> childAllele2; 

					allCountedData.addCaseMotherDuoPhased(*snpNo, reverseAlleles, motherAllele1, motherAllele2, childAllele1, childAllele2); //will divide by total later to get ave
				};

				//loop thro' case/father duos
				for(unsigned int fatherDuoNo = 1; fatherDuoNo <= noFoundFatherDuos; ++fatherDuoNo)
				{
					readHapsFile >> fatherAllele1 >> fatherAllele2 >> childAllele1 >> childAllele2; 

					allCountedData.addCaseFatherDuoPhased(*snpNo, reverseAlleles, fatherAllele1, fatherAllele2, childAllele1, childAllele2); //will divide by total later to get ave
				};

				snpNamePrev = snpName;
			};

			readHapsFile.close();
		};	

		//add to uncounted groups
		pedigreesNotCounted += (noOfSnps - includedSNPsInShapeIT.size())*noFoundTrios + 
							   (noOfSnps - includedSNPsInShapeIT.size())*noFoundMotherDuos +
							   (noOfSnps - includedSNPsInShapeIT.size())*noFoundFatherDuos;

		if(adjustedCounts) allCountedData.adjustPhasedDuoData(countAlleles, noOfSnps);

		//remove temporary files
		string rmFiles[9] = {"rm tempOutSHAPEIT"+endName+".hgraph", "rm tempOutSHAPEIT"+endName+"-hgraph.ind.me",
			"rm tempOutSHAPEIT"+endName+"-hgraph.ind.mm", "rm tempOutSHAPEIT"+endName+"-hgraph.log", "rm tempOutSHAPEIT"+endName+"-hgraph.snp.me", "rm tempOutSHAPEIT"+endName+"-hgraph.snp.mm",
			"rm tempOutSHAPEIT"+endName+".haps", "rm tempOutSHAPEIT"+endName+".sample", "rm tempOutSHAPEIT"+endName+"-max.log"};
	
		for(unsigned int rm = 0; rm <= 8; ++rm) system(rmFiles[rm].c_str());
	
	}
	else 
	{
		out("\n=============================================================================================\n");
		out("WARNING: no case/parent trio, case/mother duo or case/father duo data to phase with SHAPEIT2!\n");
		if(missThres < 0.25) out("Is the maximum missing data threshold set too low?\n");
		out("=============================================================================================\n");
	}; //end of skipping if no data to phase

	//remove temporary files
	string rmFiles2[3] = {"rm tempForSHAPEIT"+endName+".bed", "rm tempForSHAPEIT"+endName+".bim", "rm tempForSHAPEIT"+endName+".fam"};
	
	for(unsigned int rm = 0; rm <= 2; ++rm) system(rmFiles2[rm].c_str());

	
	//reopen binary file for later, doing one SNP at a time	
	readBinaryGenotypeData.open(fileName.c_str(), ios::binary);
	if(!readBinaryGenotypeData.is_open())
	{
		outErr("\nCannot read binary genotype file: "); outErr(fileName); outErr("!\n");
		exit(1);
	};

	readBinaryGenotypeData.read(buffer, 3); //magic numbers
	
};

//! Outputs lines for temp .fam file for haplotype estimatation in SHAPEIT for tios.
void ProcessData::outputPooHapTrioFamFile(const unsigned int & pedId, Pedigree * pedigree, const unsigned int & childId, ofstream & outFamPooFile)
{
	Subject * child = pedigree->getSubject(childId);
	Subject * father = pedigree->getSubject(child->getFatherId());
	Subject * mother = pedigree->getSubject(child->getMotherId());

	outFamPooFile << pedId << " " << father->getId() << " 0 0 1 1\n"; 
	outFamPooFile << pedId << " " << mother->getId() << " 0 0 2 1\n";
	outFamPooFile << pedId << " " << child->getId() << " " << father->getId() << " " << mother->getId() << " 2 2\n";
};

//! Outputs lines for temp .fam file for haplotype estimatation in SHAPEIT for mother case duos.
void ProcessData::outputPooHapMotherDuoFamFile(const unsigned int & pedId, Pedigree * pedigree, const unsigned int & childId, ofstream & outFamPooFile)
{
	Subject * child = pedigree->getSubject(childId);
	Subject * mother = pedigree->getSubject(child->getMotherId());
 
	outFamPooFile << pedId << " " << mother->getId() << " 0 0 2 1\n";
	outFamPooFile << pedId << " " << child->getId() << " 0 " << mother->getId() << " 2 2\n";
};

//! Outputs lines for temp .fam file for haplotype estimatation in SHAPEIT
void ProcessData::outputPooHapFatherDuoFamFile(const unsigned int & pedId, Pedigree * pedigree, const unsigned int & childId, ofstream & outFamPooFile)
{
	Subject * child = pedigree->getSubject(childId);
	Subject * father = pedigree->getSubject(child->getFatherId());

	outFamPooFile << pedId << " " << father->getId() << " 0 0 1 1\n"; 
	outFamPooFile << pedId << " " << child->getId() << " " << father->getId() << " 0 2 2\n";
};
//! Write binary data.
void writeBedGenotype(ofstream & psBedFile, unsigned int & psBitCount, int & aBit, const bool & allele1, const bool & allele2)
{

	aBit = aBit >> 1; //shift bits to the right
	if(allele1) aBit = aBit | 128; 
	aBit = aBit >> 1;  //shift bits to the right
	if(allele2) aBit = aBit | 128;
	
	psBitCount += 2;

	//write to file if byte is finished
	if(psBitCount == 8)
	{
		//write to file
		char buffer[1];
		buffer[0] = aBit;
		psBedFile.write(buffer, 1);
		psBitCount = 0;
		aBit = 0;
	};

};

//! Write last byte and resets the Bit and bit counter
void writeLastByteBeforeNextSNP(ofstream & psBedFile, unsigned int & psBitCount, int & aBit)
{
	if(psBitCount == 0) return; //last byte may already be written

	//shift right bits over
	while(psBitCount < 8)
	{
		aBit = aBit >> 1; //shift bits to the right
		psBitCount++;
	};
	
	//write to file
	char buffer[1];
	buffer[0] = aBit;
	psBedFile.write(buffer, 1);
	psBitCount = 0;
	aBit = 0;
};

//! Outputs lines for temp .bed file for haplotype estimatation in SHAPEIT
void ProcessData::outputPooHapBedFile(string & fileName, CountAlleles * countAlleles, map<unsigned int, unsigned int> & pooHapTrios, map<unsigned int, unsigned int> & pooHapMotherDuos, map<unsigned int, unsigned int> & pooHapFatherDuos, AllPedigreeData & allPedigreeData, set<unsigned int> & includedSNPs, ofstream & outBedPooFile)
{
	unsigned int psBitCount = 0;
	int aBit = 0;
	pair<bool, bool> fatherGeno; 
	pair<bool, bool> motherGeno; 
	pair<bool, bool> childGeno; 	
	Subject * father;
	Subject * mother;
	Subject * child;
	Pedigree * pedigree;
	set<unsigned int>::const_iterator is;

	//reopen binary SNP file, reading data one SNP at a time
	readBinaryGenotypeData.open(fileName.c_str(), ios::binary);
	if(!readBinaryGenotypeData.is_open())
	{
		outErr("\nCannot read binary genotype file: "); outErr(fileName); outErr("!\n");
		exit(1);
	};

	readBinaryGenotypeData.read(buffer, 3); //magic numbers
	

	//loop thro' SNPs and add the data for each of the chosen trios only
	for(unsigned int snpID = 1; snpID <= noOfSnps; ++snpID)
	{	
		//add SNP data if adding SNP data one at a time
		bitCount = 9; //start reading data from the next byte since in SNP-major mode for the next SNP

		//update genotype data for given SNP for every subject
		for(list<Subject *>::const_iterator s = orderListOfSubjects.begin(); s != orderListOfSubjects.end(); ++s)
		{
			addGenotypeDataToSubject(*s, snpID, countAlleles);
		};

		//check if allele 1 is minor
		//countAlleles->updateIsAllele1Minor();	 //already done
		
		is = includedSNPs.find(snpID);

		if(is != includedSNPs.end())
		{

			includedSNPsInShapeIT.push_back(snpID); //keep track of included SNPs

			//loop thro' chosen trios
			for(map<unsigned int, unsigned int>::const_iterator pht = pooHapTrios.begin(); pht != pooHapTrios.end(); ++pht)
			{
	
				if(pht->second != 0) //only pedigrees which found trios
				{
					pedigree = allPedigreeData.pedigrees.find(pht->first)->second;
					child = pedigree->getSubject(pht->second);
					father = pedigree->getSubject(child->getFatherId());
					mother = pedigree->getSubject(child->getMotherId());
	
					if(father->genotypeDataPresent(snpID)) fatherGeno = father->getGenotype(snpID); else fatherGeno = make_pair(true, false); //missing data code 				
					if(mother->genotypeDataPresent(snpID)) motherGeno = mother->getGenotype(snpID); else motherGeno = make_pair(true, false);
					if(child->genotypeDataPresent(snpID)) childGeno = child->getGenotype(snpID); else childGeno = make_pair(true, false);

					//write to .bed file in the same order as .fam file
					writeBedGenotype(outBedPooFile, psBitCount, aBit, fatherGeno.first, fatherGeno.second);
					writeBedGenotype(outBedPooFile, psBitCount, aBit, motherGeno.first, motherGeno.second);
					writeBedGenotype(outBedPooFile, psBitCount, aBit, childGeno.first, childGeno.second);
				};

			};

			//loop thro' chosen mother/case duos
			for(map<unsigned int, unsigned int>::const_iterator phmcd = pooHapMotherDuos.begin(); phmcd != pooHapMotherDuos.end(); ++phmcd)
			{
				if(phmcd->second != 0) //only pedigrees which found mother/case duos
				{
					pedigree = allPedigreeData.pedigrees.find(phmcd->first)->second;
					child = pedigree->getSubject(phmcd->second);					
					mother = pedigree->getSubject(child->getMotherId());
			
					if(mother->genotypeDataPresent(snpID)) motherGeno = mother->getGenotype(snpID); else motherGeno = make_pair(true, false); //missing data code 	
					if(child->genotypeDataPresent(snpID)) childGeno = child->getGenotype(snpID); else childGeno = make_pair(true, false);

					//write to .bed file in the same order as .fam file
					writeBedGenotype(outBedPooFile, psBitCount, aBit, motherGeno.first, motherGeno.second);
					writeBedGenotype(outBedPooFile, psBitCount, aBit, childGeno.first, childGeno.second);
				};

			};

			//loop thro' chosen father/case duos
			for(map<unsigned int, unsigned int>::const_iterator phfcd = pooHapFatherDuos.begin(); phfcd != pooHapFatherDuos.end(); ++phfcd)
			{
				if(phfcd->second != 0) //only pedigrees which found father/case duos
				{
					pedigree = allPedigreeData.pedigrees.find(phfcd->first)->second;
					child = pedigree->getSubject(phfcd->second);					
					father = pedigree->getSubject(child->getFatherId());
			
					if(father->genotypeDataPresent(snpID)) fatherGeno = father->getGenotype(snpID); else fatherGeno = make_pair(true, false); //missing data code 	
					if(child->genotypeDataPresent(snpID)) childGeno = child->getGenotype(snpID); else childGeno = make_pair(true, false);

					//write to .bed file in the same order as .fam file
					writeBedGenotype(outBedPooFile, psBitCount, aBit, fatherGeno.first, fatherGeno.second);
					writeBedGenotype(outBedPooFile, psBitCount, aBit, childGeno.first, childGeno.second);
				};

			};	

			//start new byte as starting a new SNP
			writeLastByteBeforeNextSNP(outBedPooFile, psBitCount, aBit);
		};
	};

	readBinaryGenotypeData.close();
};

