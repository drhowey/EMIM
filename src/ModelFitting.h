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


/*! \file ModelFitting.h
    \brief This file contains model fitting methods.
    
*/

#ifndef __MODELFITTING
#define __MODELFITTING

#include <map>

using namespace std;

#include "Model.h"

map<unsigned int, double> getSolnMatrixEqun(map<unsigned int, map<unsigned int, double> > & matrix, map<unsigned int, double> & vect);
	
//class to find a good fit of the parameters to the data 
class FindFit
{
private:

	Model * modelToFit; //model	
	
public:

	FindFit(Model * m) : modelToFit(m) { };

	~FindFit()
	{
		
	};

	double evaluate() const {return modelToFit->negLogLikelihood();};
	bool newtonsMethod(double & eval, double accuracy = 1e-6, unsigned int maxIterations = 50) const;
	//void displayVector2(const map<unsigned int, double> & vect) const;
	
private:

	//Newton's method for optimization methods
	void getNextPoint(map<unsigned int, double> & point) const;
	map<unsigned int, double> getPointChange() const;	
	map<unsigned int, double> getSumParameters(map<unsigned int, map<unsigned int, double> > & simplex) const;

};


#endif
