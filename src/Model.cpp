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



/*! \file Model.cpp
    \brief This file contains the source of models used for logistic regression.
    
*/

#include <map>
#include <iostream>
#include <math.h>

using namespace std; // initiates the "std" or "standard" namespace

#include "Model.h"
#include "Data.h"
#include "ModelFitting.h"

//! Sets the values of all of the parameters.
void Model::setNewParameters(map<unsigned int, double> & paras)
{
	for(map<unsigned int, double>::const_iterator p = paras.begin(); p != paras.end(); ++p)
	{
		setParameter(p->first, p->second);
	};
};

//
////! Sets the value of a given parameter for the model.
//void ModelDuoAdjustOLD::setParameter(const unsigned int & no, const double & value)
//{
//	parameters[no] = value;
//
//	switch(no)
//	{
//		case 1:
//			beta0 = value; return;
//		case 2:
//			beta1 = value; return;
//		case 3:
//			beta2 = value; return;
//		case 4:
//			beta3 = value; return;	
//	};
//};	
//
////! Finds neg for adjusting the 4a and 4b counts of duos.
//double ModelDuoAdjustOLD::negLogLikelihood()
//{
//	unsigned int noSNPs = mafs.size();
//	map<unsigned int, unsigned int>::const_iterator noFa = nosFromFather.begin(); //cell 4a
//	map<unsigned int, unsigned int>::const_iterator noMo = nosFromMother.begin(); //cell 4b
//	map<unsigned int, double>::const_iterator maf = mafs.begin(); //mafs
//
//	double negLogLik = 0;
//
//	double noDuos = (double)(noOfDuos);
//	double pq, pFa, pMo, mafSq, meanFa, meanMo, varMo, varFa, adj, adjFaDiff, adjMoDiff;
//
//	for(unsigned int i = 1; i <= noSNPs; ++i)
//	{
//		//ensure correct mafs are used
//		while(noFa->first != maf->first && noFa != nosFromFather.end()) ++noFa; 
//		while(noMo->first != maf->first && noMo != nosFromMother.end()) ++noMo; 
//
//		if(noFa != nosFromFather.end() && noMo != nosFromFather.end()) //maybe the last SNPs are missing
//		{
//
//			pq = (maf->second)*(1-maf->second);
//			if(mother)
//			{
//				pFa = (maf->second)*pq; // p^2 (1-p)
//				pMo = (1 - maf->second)*pq; // p (1-p)^2
//			}
//			else
//			{
//				pFa = (1 - maf->second)*pq; // p (1-p)^2 
//				pMo = (maf->second)*pq; // p^2 (1-p)
//			};
//			meanFa = noDuos*pFa; //N p^2 (1-p) or N p (1-p)^2 
//			meanMo = noDuos*pMo; 
//			varFa = noDuos*pFa*(1-pFa); 
//			varMo = noDuos*pMo*(1-pMo); 
//			mafSq = (maf->second)*(maf->second);
//			adj = noDuos*(beta0 + beta1*(maf->second) + beta2*mafSq + beta3*mafSq*(maf->second));
//			adjFaDiff = noFa->second - adj - meanFa;
//			adjMoDiff = noMo->second + adj - meanMo;
//
//			if(varFa > 0 && varMo > 0) negLogLik += 0.5*((adjFaDiff*adjFaDiff)/varFa + (adjMoDiff*adjMoDiff)/varMo);
//
//		};
//
//		++maf;
//		++noFa;
//		++noMo;
//	};
//
//	
//	return negLogLik;
//};
//
////! Get gradient vector for 2 parameter model.
//void ModelDuoAdjustOLD::getGradientVector(map<unsigned int, double> & gradientVector)  const
//{
//	double ans[4] = {0, 0, 0, 0};
//
//	unsigned int noSNPs = mafs.size();
//	map<unsigned int, unsigned int>::const_iterator noFa = nosFromFather.begin(); //cell 4a
//	map<unsigned int, unsigned int>::const_iterator noMo = nosFromMother.begin(); //cell 4b
//	map<unsigned int, double>::const_iterator maf = mafs.begin(); //mafs
//
//	double noDuos = (double)(noOfDuos);
//	double pq, pFa, pMo, mafSq, meanFa, meanMo, varMo, varFa, adj, adjFaDiff, adjMoDiff;
//	double bit;
//
//	for(unsigned int i = 1; i <= noSNPs; ++i)
//	{	
//		while(noFa->first != maf->first && noFa != nosFromFather.end()) ++noFa; 
//		while(noMo->first != maf->first && noMo != nosFromMother.end()) ++noMo; 
//
//		if(noFa != nosFromFather.end() && noMo != nosFromFather.end()) //maybe the last SNPs are missing
//		{
//
//			pq = (maf->second)*(1-maf->second);
//			if(mother)
//			{
//				pFa = (maf->second)*pq; // p^2 (1-p)
//				pMo = (1 - maf->second)*pq; // p (1-p)^2
//			}
//			else
//			{
//				pFa = (1 - maf->second)*pq; // p (1-p)^2 
//				pMo = (maf->second)*pq; // p^2 (1-p)
//			};
//
//			meanFa = noDuos*pFa; //N p^2 (1-p) or N p (1-p)^2 
//			meanMo = noDuos*pMo; 
//			varFa = noDuos*pFa*(1-pFa); 
//			varMo = noDuos*pMo*(1-pMo); 
//			mafSq = (maf->second)*(maf->second);
//			adj = noDuos*(beta0 + beta1*(maf->second) + beta2*mafSq + beta3*mafSq*(maf->second));
//			adjFaDiff = noFa->second - adj - meanFa;
//			adjMoDiff = noMo->second + adj - meanMo;
//
//			if(varFa > 0 && varMo > 0)
//			{
//				bit = noDuos*(-adjFaDiff/varFa + adjMoDiff/varMo);
//				ans[0] += bit;
//				ans[1] += bit*(maf->second);
//				ans[2] += bit*mafSq;
//				ans[3] += bit*mafSq*(maf->second);
//			};
//
//		};
//
//		++maf;
//		++noFa;
//		++noMo;
//	};
//
//	gradientVector[1] = ans[0];	
//	gradientVector[2] = ans[1];
//	gradientVector[3] = ans[2];	
//	gradientVector[4] = ans[3];	
//};
//
////! Returns the 2nd derivative w.r.t. chosen parameters of the negative log likelihood.
//void ModelDuoAdjustOLD::getHessianMatrix(map<unsigned int, map<unsigned int, double> > & hessianMatrix) const
//{
//	
//	double ans[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}; //col, row
//	
//	unsigned int noSNPs = mafs.size();
//	map<unsigned int, unsigned int>::const_iterator noFa = nosFromFather.begin(); //cell 4a
//	map<unsigned int, unsigned int>::const_iterator noMo = nosFromMother.begin(); //cell 4b
//	map<unsigned int, double>::const_iterator maf = mafs.begin(); //mafs
//
//	double noDuos = (double)(noOfDuos);
//	double pq, pFa, pMo, mafSq, maf3, maf4, maf5, varMo, varFa;
//	double bit;
//
//	for(unsigned int i = 1; i <= noSNPs; ++i)
//	{	
//		while(noFa->first != maf->first && noFa != nosFromFather.end()) ++noFa; 
//		while(noMo->first != maf->first && noMo != nosFromMother.end()) ++noMo; 
//
//		if(noFa != nosFromFather.end() && noMo != nosFromFather.end()) //maybe the last SNPs are missing
//		{
//
//			pq = (maf->second)*(1-maf->second);
//			if(mother)
//			{
//				pFa = (maf->second)*pq; // p^2 (1-p)
//				pMo = (1 - maf->second)*pq; // p (1-p)^2
//			}
//			else
//			{
//				pFa = (1 - maf->second)*pq; // p (1-p)^2 
//				pMo = (maf->second)*pq; // p^2 (1-p)
//			};
//			varFa = noDuos*pFa*(1-pFa); 
//			varMo = noDuos*pMo*(1-pMo); 
//			mafSq = (maf->second)*(maf->second);
//			maf3 = mafSq*(maf->second);
//			maf4 = maf3*(maf->second);
//			maf5 = maf4*(maf->second);
//		
//			if(varFa > 0 && varMo > 0)
//			{
//				bit = noDuos*noDuos*(1/varFa + 1/varMo);
//				ans[0][0] += bit;		
//				ans[0][1] += bit*(maf->second);	
//				ans[0][2] += bit*mafSq;
//				ans[0][3] += bit*maf3;
//
//				ans[1][1] += bit*mafSq;
//				ans[1][2] += bit*maf3;
//				ans[1][3] += bit*maf4;
//
//				ans[2][2] += bit*maf4;
//				ans[2][3] += bit*maf5;
//
//				ans[3][3] += bit*maf5*(maf->second);
//			};
//
//		};
//
//		++maf;
//		++noFa;
//		++noMo;
//	};
//		
//	//The Hessian matrix is symetric so only calculate half and then copy	
//	ans[1][0] = ans[0][1];	
//	ans[2][0] = ans[0][2];	
//	ans[2][1] = ans[1][2];	
//	ans[3][0] = ans[0][3];	
//	ans[3][1] = ans[1][3];	
//	ans[3][2] = ans[2][3];	
//	
//	//setup the matrix with calculated values
//	map<unsigned int, double> aCol;	
//	for(unsigned int col = 0; col < 4; ++col)
//	{		
//		for(unsigned int row = 0; row < 4; ++row)
//		{
//			aCol[row + 1] = ans[col][row];	
//		};
//		hessianMatrix[col + 1] = aCol;
//	};
//
//};
//

//! Sets the value of a given parameter for the model.
void ModelDuoAdjust::setParameter(const unsigned int & no, const double & value)
{
	parameters[no] = value;

	switch(no)
	{
		case 1:
			beta0 = value; return;
		case 2:
			beta1 = value; return;
		case 3:
			beta2 = value; return;
		case 4:
			beta3 = value; return;	
	};
};	

//! Finds neg for adjusting the 4a and 4b counts of duos.
double ModelDuoAdjust::negLogLikelihood()
{
	unsigned int noSNPs = mafs.size();
	map<unsigned int, unsigned int>::const_iterator noFa = nosFromFather.begin(); //cell 4a
	map<unsigned int, unsigned int>::const_iterator noMo = nosFromMother.begin(); //cell 4b
	map<unsigned int, double>::const_iterator maf = mafs.begin(); //mafs

	double negLogLik = 0;

	double noDuos = (double)(noOfDuos);
	double pq, pFa, pMo, mafSq, meanFa, meanMo, varMo, varFa, covarFaMo, adj, adjFaDiff, adjMoDiff;

	for(unsigned int i = 1; i <= noSNPs; ++i)
	{
		//ensure correct mafs are used
		while(noFa->first != maf->first && noFa != nosFromFather.end()) ++noFa; 
		while(noMo->first != maf->first && noMo != nosFromMother.end()) ++noMo; 

		if(noFa != nosFromFather.end() && noMo != nosFromFather.end()) //maybe the last SNPs are missing
		{

			pq = (maf->second)*(1-maf->second);
			if(mother)
			{
				pFa = (maf->second)*pq; // p^2 (1-p)
				pMo = (1 - maf->second)*pq; // p (1-p)^2
			}
			else
			{
				pFa = (1 - maf->second)*pq; // p (1-p)^2 
				pMo = (maf->second)*pq; // p^2 (1-p)
			};

			meanFa = noDuos*pFa; //N p^2 (1-p) or N p (1-p)^2 
			meanMo = noDuos*pMo; 
			mafSq = (maf->second)*(maf->second);
			adj = noDuos*(beta0 + beta1*(maf->second) + beta2*mafSq + beta3*mafSq*(maf->second));
			adjFaDiff = noFa->second - adj - meanFa;
			adjMoDiff = noMo->second + adj - meanMo;
			varFa = noDuos*pFa*(1-pFa); 
			varMo = noDuos*pMo*(1-pMo);
			covarFaMo = -noDuos*pFa*pMo;

			if(pq > 0) negLogLik += 0.5*(adjFaDiff*adjFaDiff*varMo - 2*adjFaDiff*adjMoDiff*covarFaMo + adjMoDiff*adjMoDiff*varFa)/(varFa*varMo - covarFaMo*covarFaMo);
		};

		++maf;
		++noFa;
		++noMo;
	};

	return negLogLik;
};

//! Get gradient vector for 2 parameter model.
void ModelDuoAdjust::getGradientVector(map<unsigned int, double> & gradientVector)  const
{
	double ans[4] = {0, 0, 0, 0};

	unsigned int noSNPs = mafs.size();
	map<unsigned int, unsigned int>::const_iterator noFa = nosFromFather.begin(); //cell 4a
	map<unsigned int, unsigned int>::const_iterator noMo = nosFromMother.begin(); //cell 4b
	map<unsigned int, double>::const_iterator maf = mafs.begin(); //mafs

	double noDuos = (double)(noOfDuos);
	double pq, pFa, pMo, mafSq, meanFa, meanMo, varMo, varFa, covarFaMo, adj, adjFaDiff, adjMoDiff;
	double bit;

	for(unsigned int i = 1; i <= noSNPs; ++i)
	{	
		while(noFa->first != maf->first && noFa != nosFromFather.end()) ++noFa; 
		while(noMo->first != maf->first && noMo != nosFromMother.end()) ++noMo; 

		if(noFa != nosFromFather.end() && noMo != nosFromFather.end()) //maybe the last SNPs are missing
		{

			pq = (maf->second)*(1-maf->second);
			if(mother)
			{
				pFa = (maf->second)*pq; // p^2 (1-p)
				pMo = (1 - maf->second)*pq; // p (1-p)^2
			}
			else
			{
				pFa = (1 - maf->second)*pq; // p (1-p)^2 
				pMo = (maf->second)*pq; // p^2 (1-p)
			};

			meanFa = noDuos*pFa; //N p^2 (1-p) or N p (1-p)^2 
			meanMo = noDuos*pMo; 
			varFa = noDuos*pFa*(1-pFa); 
			varMo = noDuos*pMo*(1-pMo); 
			mafSq = (maf->second)*(maf->second);
			adj = noDuos*(beta0 + beta1*(maf->second) + beta2*mafSq + beta3*mafSq*(maf->second));
			adjFaDiff = noFa->second - adj - meanFa;
			adjMoDiff = noMo->second + adj - meanMo;
			covarFaMo = -noDuos*pFa*pMo;
		
		    if(pq > 0)
			{
				bit = noDuos*(-adjFaDiff*(varMo + covarFaMo) + adjMoDiff*(varFa + covarFaMo))/(varFa*varMo - covarFaMo*covarFaMo);
				ans[0] += bit;
				ans[1] += bit*(maf->second);
				ans[2] += bit*mafSq;
				ans[3] += bit*mafSq*(maf->second);
			};

		};

		++maf;
		++noFa;
		++noMo;
	};

	gradientVector[1] = ans[0];	
	gradientVector[2] = ans[1];
	gradientVector[3] = ans[2];	
	gradientVector[4] = ans[3];	

};

//! Returns the 2nd derivative w.r.t. chosen parameters of the negative log likelihood.
void ModelDuoAdjust::getHessianMatrix(map<unsigned int, map<unsigned int, double> > & hessianMatrix) const
{
	
	double ans[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}; //col, row
	
	unsigned int noSNPs = mafs.size();
	map<unsigned int, unsigned int>::const_iterator noFa = nosFromFather.begin(); //cell 4a
	map<unsigned int, unsigned int>::const_iterator noMo = nosFromMother.begin(); //cell 4b
	map<unsigned int, double>::const_iterator maf = mafs.begin(); //mafs

	double noDuos = (double)(noOfDuos);
	double pq, pFa, pMo, mafSq, maf3, maf4, maf5, varMo, varFa, covarFaMo;
	double bit;

	for(unsigned int i = 1; i <= noSNPs; ++i)
	{	
		while(noFa->first != maf->first && noFa != nosFromFather.end()) ++noFa; 
		while(noMo->first != maf->first && noMo != nosFromMother.end()) ++noMo; 

		if(noFa != nosFromFather.end() && noMo != nosFromFather.end()) //maybe the last SNPs are missing
		{

			pq = (maf->second)*(1-maf->second);
			if(mother)
			{
				pFa = (maf->second)*pq; // p^2 (1-p)
				pMo = (1 - maf->second)*pq; // p (1-p)^2
			}
			else
			{
				pFa = (1 - maf->second)*pq; // p (1-p)^2 
				pMo = (maf->second)*pq; // p^2 (1-p)
			};

			varFa = noDuos*pFa*(1-pFa); 
			varMo = noDuos*pMo*(1-pMo); 
			covarFaMo = -noDuos*pFa*pMo;
			mafSq = (maf->second)*(maf->second);
			maf3 = mafSq*(maf->second);
			maf4 = maf3*(maf->second);
			maf5 = maf4*(maf->second);
	
			if(pq > 0)
			{	
				bit = noDuos*noDuos*(varFa + 2*covarFaMo + varMo)/(varFa*varMo - covarFaMo*covarFaMo);

				ans[0][0] += bit;		
				ans[0][1] += bit*(maf->second);	
				ans[0][2] += bit*mafSq;
				ans[0][3] += bit*maf3;

				ans[1][1] += bit*mafSq;
				ans[1][2] += bit*maf3;
				ans[1][3] += bit*maf4;

				ans[2][2] += bit*maf4;
				ans[2][3] += bit*maf5;

				ans[3][3] += bit*maf5*(maf->second);
			};

		};

		++maf;
		++noFa;
		++noMo;
	};
		
	//The Hessian matrix is symetric so only calculate half and then copy	
	ans[1][0] = ans[0][1];	
	ans[2][0] = ans[0][2];	
	ans[2][1] = ans[1][2];	
	ans[3][0] = ans[0][3];	
	ans[3][1] = ans[1][3];	
	ans[3][2] = ans[2][3];	
	
	//setup the matrix with calculated values
	map<unsigned int, double> aCol;	
	for(unsigned int col = 0; col < 4; ++col)
	{		
		for(unsigned int row = 0; row < 4; ++row)
		{
			aCol[row + 1] = ans[col][row];	
		};
		hessianMatrix[col + 1] = aCol;
	};

};

