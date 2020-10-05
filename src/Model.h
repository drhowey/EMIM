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


/*! \file Model.h
    \brief This file contains the models used for logistic regression.
    
*/

#ifndef __MODEL
#define __MODEL


#include <map>
#include <list>

using namespace std;

struct JointGenotypeCounts;
struct SNPData;
struct QuantitiveTraits;
struct CovariateData;

#include "Data.h"

//! General Class for Models.
class Model
{
protected:

	map<unsigned int, double> parameters;
	
	
public:

	Model() : parameters() {};
	//setup parameter values and initValues of variables, to be over written in subclass of each model
	//where the variables created and added to the vector for the variables
	Model(map<unsigned int, double> paras) : parameters(paras) {};
	
	virtual ~Model() {};

	double getParameter(const unsigned int & no) const {return parameters.find(no)->second;};
	map<unsigned int, double> getParameters() const {return parameters;};
	void setNewParameters(map<unsigned int, double> & paras);

	
	virtual void setParameter(const unsigned int & no, const double & value) {parameters[no] = value;};	
	virtual double negLogLikelihood() {return 0;};
	virtual void getGradientVector(map<unsigned int, double> & gradientVector) const {};
	virtual void getHessianMatrix(map<unsigned int, map<unsigned int, double> > & hessianMatrix) const {};

};

//! General Class for Models with 2 SNPs.
class ModelDuoAdjust : public Model
{
private:
	
	double beta0, beta1, beta2, beta3;
	bool mother;

	map<unsigned int, unsigned int> nosFromFather; //snpID, cell 4a count, minor alleles from father *not* risk alleles
	map<unsigned int, unsigned int> nosFromMother; //snpID, cell 4b count, minor alleles from mother *not* risk alleles
	map<unsigned int, double> mafs; //mafs
	unsigned int noOfDuos; //the same for all SNPs as given to SHAPEIT2 together

public:

	ModelDuoAdjust(const bool & m) : mother(m) {};

	//setup parameter values and initValues of variables, to be over written in subclass of each model
	//where the variables created and added to the vector for the variables
		
	~ModelDuoAdjust() {};

	void setParameter(const unsigned int & no, const double & value);
	double negLogLikelihood();
	void getGradientVector(map<unsigned int, double> & gradientVector) const;
	void getHessianMatrix(map<unsigned int, map<unsigned int, double> > & hessianMatrix) const;

	void clearData() {nosFromFather.clear(); nosFromMother.clear();};
	void setNoOfDuos(const unsigned int & nd) {noOfDuos = nd;};
	void addFromFatherData(const unsigned int & snpID, const unsigned int & no) {nosFromFather[snpID] = no;};
	void addFromMotherData(const unsigned int & snpID, const unsigned int & no) {nosFromMother[snpID] = no;};
	void addMAFData(map<unsigned int, double> & m) {mafs = m;};
};


////! General Class for Models with 2 SNPs.
//class ModelDuoAdjustOLD : public Model
//{
//private:
//	
//	double beta0, beta1, beta2, beta3;
//	bool mother;
//
//	map<unsigned int, unsigned int> nosFromFather; //snpID, cell 4a count, minor alleles from father *not* risk alleles
//	map<unsigned int, unsigned int> nosFromMother; //snpID, cell 4b count, minor alleles from mother *not* risk alleles
//	map<unsigned int, double> mafs; //mafs
//	unsigned int noOfDuos; //the same for all SNPs as given to SHAPEIT2 together
//
//public:
//
//	ModelDuoAdjustOLD(const bool & m) : mother(m) {};
//
//	//setup parameter values and initValues of variables, to be over written in subclass of each model
//	//where the variables created and added to the vector for the variables
//		
//	~ModelDuoAdjustOLD() {};
//
//	void setParameter(const unsigned int & no, const double & value);
//	double negLogLikelihood();
//	void getGradientVector(map<unsigned int, double> & gradientVector) const;
//	void getHessianMatrix(map<unsigned int, map<unsigned int, double> > & hessianMatrix) const;
//
//	void clearData() {nosFromFather.clear(); nosFromMother.clear();};
//	void setNoOfDuos(const unsigned int & nd) {noOfDuos = nd;};
//	void addFromFatherData(const unsigned int & snpID, const unsigned int & no) {nosFromFather[snpID] = no;};
//	void addFromMotherData(const unsigned int & snpID, const unsigned int & no) {nosFromMother[snpID] = no;};
//	void addMAFData(map<unsigned int, double> & m) {mafs = m;};
//};

#endif
