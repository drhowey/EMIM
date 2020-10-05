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


#ifndef __MAIN
#define __MAIN

extern bool outputToScreen; 

extern ofstream logFile; 

template<typename T>
//! Outputs message to screen and log file
void out(const T & text)
{	
	if(outputToScreen) cout << text;
	logFile << text;
};

template<typename T>
//! Outputs error message to screen and log file
void outErr(const T & text)
{
	cerr << text;
	logFile << text;
};


//! This defines a compiler variable to whether gzipped files are handled by PREMIM.
//! Comment out this line if gzip files are not required. If this variable is set then
//! the zlib and gzstream libraries are required.
//#define USING_GZIP

#endif

/*! \mainpage PREMIM Source Code Documenation
 *
 * \section intro_sec Introduction
 *
 * This documentation is automatically produced from the source code comments.
 * Hopefully it will be useful in providing an overview of how PREMIM works. 
 *  
 */

