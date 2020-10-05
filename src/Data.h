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


/*! \file Data.h
    \brief This file contains basic classes for storing genotype infomation.
    
*/

#ifndef __DATA
#define __DATA

#include <string>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <ostream>

#include "Model.h"

using namespace std; // initiates the "std" or "standard" namespace

//! Changes the marker output file if necessary.
void updateMarkersFile(ofstream & markersFile, const unsigned int & snpId, const unsigned int & splitSNPOutput, const string & outputDirectory, const string & endName);

//! Adds to counts.
void addToChildrenCountForTrioOrDuo(map<unsigned int, map<unsigned int, unsigned int> > & childrenCount, const unsigned int & pedId, const unsigned int & childId);

//! Maps the names of subjects and pedigrees given in file to an ordinal number.
class MapIds
{
private:
	map<string, unsigned int> idMap; //!< original name in file, counted id

public:

	MapIds() : idMap() {};

	//! Create a new ID for an item.
	unsigned int addItem(const string & name)
	{
		unsigned int newId = (unsigned)idMap.size() + 1;
		idMap[name] = newId;
		return newId; 
	};

	//! Look up ID from a name given from file.
	unsigned int getId(const string & name);
	//! Look up name from a ID
	string getName(const unsigned int & id);
};

//! Maps the allele names given in file to either a 0 or 1.

//! In the .ped pedigree file the SNP names may be given by a letter, e.g. A and G,
//! these letters are mapped to values 0 and 1 in the code which in turn map to
//! allele "1" and allele "2" respectively.
class MapSnpAlleleIds
{
private:
	map<unsigned int, map<string, bool> > mapSnpAlleles; //SNP ID no., allele name,  allele id - 0 or 1

public:

	MapSnpAlleleIds() {};

	~MapSnpAlleleIds() {};

	bool getSnpAlleleId(const unsigned int & snpId, string & alleleName); 
	map<string, bool> getSnpAlleleNames(const unsigned int & snpId) const; 
};

//! The genotype data of one subject.
class GenotypeData
{
private:
	map<unsigned int, pair<bool, bool> > data;  //snp ID, (allele 1, allele 2), 0 = allele value 1, 1 = allele value 2

public:

	GenotypeData() : data()
	{
	  
	};

	~GenotypeData()
	{
		data.clear();
	};

	void addGenotype(const unsigned int & snpId, const bool & allele1, const bool & allele2);
	pair<bool, bool> getGenotype(const unsigned int & snpId) const;
	bool genotypeDataPresent(const unsigned int & snpId) const;
};

//! A subject with its attributes including parents.
class Subject
{
private:
	unsigned int id;
	unsigned int fatherId;
	unsigned int motherId;
	unsigned int sex; //1 = male, 2 = female, 0 = unknown
	bool affected;
	bool removed;
	
public:

	Subject(unsigned int i, unsigned int fi, unsigned int mi, unsigned int s, bool a): id(i), fatherId(fi), motherId(mi), sex(s), affected(a), removed(false) {};

	virtual ~Subject() {};
	
	virtual void addGenotype(const unsigned int & snpId, const bool & all1, const bool & all2) {};
	virtual pair<bool, bool> getGenotype(const unsigned int & snpId) const {return make_pair(true, false);};
	virtual bool genotypeDataPresent(const unsigned int & snpId) const {return false;};
	
	bool getAffected() const {return affected;};
	unsigned int getFatherId() const {return fatherId;};
	unsigned int getMotherId() const {return motherId;};
	unsigned int getId() const {return id;};
	void removeFromPedigree() {removed = true;};
	void restoreToPedigree() {removed = false;};
	bool isNotRemoved() {return !removed;};
	void outputDetails();
};

//! A subject including all of its genotype data.
class SubjectMultiSnp : public Subject
{
private:
	GenotypeData * genotypeData;

public:

	SubjectMultiSnp(unsigned int i, unsigned int fi, unsigned int mi, unsigned int s, bool a): Subject(i, fi, mi, s, a)
	{
		genotypeData = new GenotypeData();
	};

	~SubjectMultiSnp()
	{
		delete genotypeData;		
	};

	void addGenotype(const unsigned int & snpId, const bool & all1, const bool & all2) {genotypeData->addGenotype(snpId, all1, all2);};
	pair<bool, bool> getGenotype(const unsigned int & snpId) const {return genotypeData->getGenotype(snpId);};
	bool genotypeDataPresent(const unsigned int & snpId) const {return genotypeData->genotypeDataPresent(snpId);};
	
};

//! A subject with genotype data for only one SNP.
class SubjectOneSnp : public Subject
{
private:
	bool allele1, allele2;

public:

	SubjectOneSnp(unsigned int i, unsigned int fi, unsigned int mi, unsigned int s, bool a): Subject(i, fi, mi, s, a) {};
	
	~SubjectOneSnp() {};

	void addGenotype(const unsigned int & snpId, const bool & all1, const bool & all2) {allele1 = all1; allele2 = all2;};
	pair<bool, bool> getGenotype(const unsigned int & snpId) const
	{		
		return make_pair(allele1, allele2);
	};
	bool genotypeDataPresent(const unsigned int & snpId) const {return (!allele1 || allele2);}; //genotype 1 / 0 denotes missing genotype	
};

//! A simple grouping of a father, mother and child.
struct Trio
{
	Subject * father;
	Subject * mother;
	Subject * child;

	Trio() : father(0), mother(0), child(0) {};
	~Trio() {};
};

//! Contains a list of all subjects belonging to a pedigree.
class Pedigree
{
private:
	unsigned int pedId;
	map<unsigned int, Subject *> subjects;//subject Id, pointer to subject object
	bool pooHapProcessed;

public:

	Pedigree(unsigned int & p) : pedId(p), subjects(), pooHapProcessed(false) {};

	~Pedigree()
	{
		for(map<unsigned int, Subject *>::iterator s = subjects.begin(); s != subjects.end(); ++s)
		{
			delete s->second;
		};
	};

	void addSubject(unsigned int & subjectId, Subject * subject)
	{
		subjects[subjectId] = subject;
	};

	unsigned int getPedId() {return pedId;};
	Subject * getSubject(const unsigned int & id) const;
	Subject * getSubject(const unsigned int & id, bool & exists) const;
	Trio findCaseParentTrio(const unsigned int & snpId, bool & found);
	pair<Subject *, Subject *> findCaseMotherDuo(const unsigned int & snpId, bool & found);
	pair<Subject *, Subject *> findCaseFatherDuo(const unsigned int & snpId, bool & found);
	Subject * findCase(const unsigned int & snpId, bool & found);
	pair<Subject *, Subject *> findParentsOfCase(const unsigned int & snpId, bool & found);
	Subject * findMotherOfCase(const unsigned int & snpId, bool & found);
	Subject * findFatherOfCase(const unsigned int & snpId, bool & found);
	Trio findParentsOfControl(const unsigned int & snpId, bool & found);
	pair<Subject *, Subject *> findControlMotherDuo(const unsigned int & snpId, bool & found);
	pair<Subject *, Subject *> findControlFatherDuo(const unsigned int & snpId, bool & found);
	Subject * findControl(const unsigned int & snpId, bool & found);

	void findCaseParentTrio(Trio & aTrio, const unsigned int & snpId, bool & found, Subject * subject);
	void findCaseMotherDuo(pair<Subject *, Subject *> & aCaseMotherDuo, const unsigned int & snpId, bool & found, Subject * subject);
	void findCaseFatherDuo(pair<Subject *, Subject *> & aCaseFatherDuo, const unsigned int & snpId, bool & found, Subject * subject);
	void findCase(Subject *& aCase, const unsigned int & snpId, bool & found, Subject * subject);
	void findParentsOfCase(pair<Subject *, Subject *> & parents, const unsigned int & snpId, bool & found, Subject * subject);
	void findMotherOfCase(Subject *& mother, const unsigned int & snpId, bool & found, Subject * subject);
	void findFatherOfCase(Subject *& father, const unsigned int & snpId, bool & found, Subject * subject);
	void findParentsOfControl(Trio & aTrio, const unsigned int & snpId, bool & found, Subject * subject);
	void findControlMotherDuo(pair<Subject *, Subject *> & aControlMotherDuo, const unsigned int & snpId, bool & found, Subject * subject);
	void findControlFatherDuo(pair<Subject *, Subject *> & aControlFatherDuo, const unsigned int & snpId, bool & found, Subject * subject);
	void findControl(Subject *& control, const unsigned int & snpId, bool & found, Subject * subject);

	bool findPooHapCaseParentTrios(const unsigned int & snpId, map<unsigned int, map<unsigned int, unsigned int> > & childrenCountForTrio);
	bool findPooHapCaseMotherDuos(const unsigned int & snpId, map<unsigned int, map<unsigned int, unsigned int> > & childrenCountForMotherDuo);
	bool findPooHapCaseFatherDuos(const unsigned int & snpId, map<unsigned int, map<unsigned int, unsigned int> > & childrenCountForFatherDuo);

	unsigned int getNumberOfSubjects() const {return (unsigned)subjects.size();};
	void restoreAllSubjectsToPedigree();
	void setPooHapProcessed(const bool & b) {pooHapProcessed = b;};
	bool isPooHapProcessed() const {return pooHapProcessed;};
}; 

//! Allele counts for SNP "1" for affected and unaffecteds
struct AlleleCounts
{
	unsigned int unaffectedAllele1;
	unsigned int affectedAllele1;
	unsigned int totalSubjectsUnaffected; //this does not include unknown genotypes, that is why it is counted per SNP
	unsigned int totalSubjectsAffected; //this does not include unknown genotypes

	AlleleCounts() : unaffectedAllele1(0), affectedAllele1(0), totalSubjectsUnaffected(0), totalSubjectsAffected(0) {};

	AlleleCounts(const unsigned int & uaa1, const unsigned int & aa1, const unsigned int & tsa, const unsigned int & tsu) : 
	unaffectedAllele1(uaa1), affectedAllele1(aa1), totalSubjectsUnaffected(tsa), totalSubjectsAffected(tsu) {};
};

//! This records whether allele "1" is minor or not
class CountAlleles
{
private:

	map<unsigned int, bool> isAllele1Minor; //snpId, is allele1 minor?
	map<unsigned int, AlleleCounts> allele1Counts; //snpId, unaffected count, affected count
	map<unsigned int, bool> wereAllelesReversed; //snpId, were Alleles reversed?

public:

	CountAlleles() : isAllele1Minor() {};

	~CountAlleles() {};

	bool isSnpAllele1Minor(const unsigned int & snpId) const;	
	void clearAllele1Counts() {allele1Counts.clear();};	
	void addToAllele1Count(const unsigned int & snpId, const unsigned int & i, const bool & affected);
	void updateIsAllele1Minor();
	void updateIsAllele1Minor(unsigned int & snpID);
	bool isSNPDataMissing(unsigned int & snpId) const;
	void addReverseAlleles(const unsigned int & snpId, bool & reverseAlleles) {wereAllelesReversed[snpId] = reverseAlleles;};
	bool getWereAllelesReversed(const unsigned int & snpId) const;
	void outputMarkersFile(ofstream & markersFile, const unsigned int & noOfSnps, const unsigned int & splitSNPOutput, const string & outputDirectory, const string & endName, bool & useMajorAlleleAsRisk, map<unsigned int, string> & riskAlleleNames, MapSnpAlleleIds & snpAlleleIds, const bool & useGivenRiskAlleles);
	void outputMarkersFileOneLine(ofstream & markersFile, const unsigned int & snpId, const unsigned int & splitSNPOutput, const string & outputDirectory, const string & endName, bool & useMajorAlleleAsRisk, map<unsigned int, string> & riskAlleleNames, MapSnpAlleleIds & snpAlleleIds, const bool & useGivenRiskAlleles);
	map<unsigned int, double> getMinorAlleleFreqs(unsigned int & noOfSNPs);
};

//! General class containing methods for counted genotype data
class CountedGenotypeData
{
private:
	map<unsigned int, map<unsigned int, unsigned int> > countedData; //SNP ID, group ID, count
	
	unsigned int totalSamplesPhased;
	unsigned int totalSnpGroups;
	unsigned int undefinedSnpGroups; //SNP-groups unplaced due to a Mendelian error
	string filename;
	unsigned int noOfGroups;
	string groupName;
	ofstream fileOutput;
	bool somethingCounted;

protected:
	map<unsigned int, map<unsigned int, unsigned int> > countedDataPhased; //SNP ID, group ID, count
	map<unsigned int, map<unsigned int, double> > countedDataPhasedAdjusted; //SNP ID, group ID, adjusted count - for adjusting duos
	bool adjustedCounts;

public:

	CountedGenotypeData(const string & fname, const unsigned int & ngrps, const string & grpNm, const string & dir, const string & en, const unsigned int & splitSNPOutput, const bool & pooHaps = false)
		: countedData(), totalSamplesPhased(0), totalSnpGroups(0), undefinedSnpGroups(0), filename(dir+fname+en),
		noOfGroups(ngrps), groupName(grpNm), somethingCounted(false)
	{
		//open files, do now in case we are writing results SNP by SNP
		if(splitSNPOutput == 0)
		{
			string snpFilename = filename+".dat";
			fileOutput.open(snpFilename.c_str());
			if(!pooHaps) fileOutput << "snp\tcellcount 1-"<<noOfGroups<<"\n";
			else fileOutput << "snp\tcellcount 1-"<<noOfGroups<<" (+2 phased, +2 totals)\n";
		};
	};

	~CountedGenotypeData()
	{
		if(fileOutput.is_open()) fileOutput.close();
	};

	void addItem(const unsigned int & snpId, unsigned int & groupId);
	void addItemPhased(const unsigned int & snpId, unsigned int & groupId);
	void outputResultsFile(const unsigned int & noOfSnps, const unsigned int & totalNoOfPedigrees, const unsigned int & splitSNPOutput);
	void outputResultFileOneSnp(const unsigned int & snpId, const unsigned int & splitSNPOutput, const bool & pooHaps = false);
	unsigned int hasData() const
	{
	 if(!countedDataPhased.empty()) return 2;
	 else if(!countedData.empty() || somethingCounted) return 1;
	 else return 0;
	};
	void writeSnpGroups(ofstream & fileOutput, const unsigned int & snpId, map<unsigned int, map<unsigned int, unsigned int> >::const_iterator & snp, map<unsigned int, map<unsigned int, unsigned int> >::const_iterator & snpPh, const bool & oneSnp);
	void writeSnpGroupsPhased(ofstream & fileOutput, const unsigned int & snpId, map<unsigned int, map<unsigned int, unsigned int> >::const_iterator & snp, map<unsigned int, map<unsigned int, unsigned int> >::const_iterator & snpPh, const bool & oneSnp);
	void outputFileSummary(const unsigned int & noOfSnps, const unsigned int & totalNoOfPedigrees);
	void setTotalSamples(unsigned int & totSams) {totalSamplesPhased = totSams;};
	void setNoOfGroups(const unsigned int & ng) {noOfGroups = ng;};
	void setAdjustedCounts(const bool & adj) {adjustedCounts = adj;};
};

//! Class for storing data on trios
class CountedTrios : public CountedGenotypeData
{
private:


public:
	CountedTrios(const string & filename, const unsigned int & ngrps, const string & grpNm, const string & dir, const string & en, const unsigned int & splitSNPOutput, bool & ph) : CountedGenotypeData(filename, ngrps, grpNm, dir, en, splitSNPOutput, ph) {};

	~CountedTrios()
	{
		
	};

	unsigned int getGroup(const unsigned int & snpId, Subject * father, Subject * mother, Subject * child);
	unsigned int getGroupPhased(const unsigned int & snpId, bool & fatherAllele1, bool & fatherAllele2, bool & motherAllele1, bool & motherAllele2, bool & childAllele1, bool & childAllele2);
	void addTrio(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * father, Subject * mother, Subject * child);
	void addTrioPhased(const unsigned int & snpId, const bool & reverseAlleleLabels, bool & fatherAllele1, bool & fatherAllele2, bool & motherAllele1, bool & motherAllele2, bool & childAllele1, bool & childAllele2);
	void changeGroupIdForReverseLabels(unsigned int & groupId);
};

//! Class for storing data on duos
class CountedDuos : public CountedGenotypeData
{
private:
	
public:
	CountedDuos(const string & filename, const unsigned int & ngrps, const string & grpNm, const string & dir, const string & en, const unsigned int & splitSNPOutput, const bool & ph = false) : CountedGenotypeData(filename, ngrps, grpNm, dir, en, splitSNPOutput, ph) {};

	~CountedDuos() 
	{
		
	};

	unsigned int getGroup(const unsigned int & snpId, Subject * parent, Subject * child);
	unsigned int getGroupPhased(const unsigned int & snpId, bool & parentAllele1, bool & parentAllele2, bool & childAllele1, bool & childAllele2, const bool & mother);
	void addDuo(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * parent, Subject * child);
	void addDuoPhased(const unsigned int & snpId, const bool & reverseAlleleLabels, bool & parentAllele1, bool & parentAllele2, bool & childAllele1, bool & childAllele2, const bool & mother);
	void changeGroupIdForReverseLabels(unsigned int & groupId);
	void adjustPhasedDuoData(CountAlleles * countAlleles, unsigned int & noOfSNPs, const bool & mother);
};

//! Class for storing data on parents
class CountedParents : public CountedGenotypeData
{
private:


public:
	CountedParents(const string & filename, const unsigned int & ngrps, const string & grpNm, const string & dir, const string & en, const unsigned int & splitSNPOutput) : CountedGenotypeData(filename, ngrps, grpNm, dir, en, splitSNPOutput) {};

	~CountedParents()
	{
		
	};

	unsigned int getGroup(const unsigned int & snpId, Subject * father, Subject * mother);
	void addParents(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * father, Subject * mother);
	void changeGroupIdForReverseLabels(unsigned int & groupId);
};

//! Class for storing data on a single subject
class CountedSingleSubject : public CountedGenotypeData
{
private:


public:
	CountedSingleSubject(const string & filename, const unsigned int & ngrps, const string & grpNm, const string & dir, const string & en, const unsigned int & splitSNPOutput) : CountedGenotypeData(filename, ngrps, grpNm, dir, en, splitSNPOutput) {};

	~CountedSingleSubject() 
	{
		
	};

	unsigned int getGroup(const unsigned int & snpId, Subject * subject);
	void addSubject(const unsigned int & snpId, const bool & reverseAlleleLabels, Subject * subject);
	void changeGroupIdForReverseLabels(unsigned int & groupId);
};



#endif
