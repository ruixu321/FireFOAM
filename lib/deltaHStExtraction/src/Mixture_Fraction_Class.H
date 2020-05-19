/***************************************************************************
 *   Copyright (C) 2018 by Van Minh LE, Alexis Marchand                    *
 *   van-minh.le@ensma.fr  amarcha1@umd.edu                                *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
/*---------------------------------------------------------------------------*\
Class
	Mixture_Fraction_Class
Description
    Mixture_Fraction_Class is a step in the whole process, which reads the 1st
 	LUT to find deltaHSt. The class is a treatment for mixture fraction Z 
	distribution
SourceFiles
    Mixture_Fraction_Class.C
\*---------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;


#ifndef Mixture_Fraction_Class_H
#define Mixture_Fraction_Class_H

#include "Variance_Mixture_Fraction_Class.H"

class Mixture_Fraction_Class
{

public:

    Mixture_Fraction_Class();           // Constructor
	void Read();                        // Reading the LUT
    void GetMeanValues(const double Z, const double Zv, const double ChiSt, const double deltaH, double &extractedDeltaHSt);

private:
    std::string         _path_library;          // path to the LUT
    std::string         _name;                  // file name, in which Z distribution information is clarified 

    int     _nZ;                                // number of mixture fraction Z
    std::vector<double>         _Z;             // Z values
    std::vector<string>         _Z_name;        // corresponding directory name to the value of Z
    Variance_Mixture_Fraction_Class         *_Zv;        // Pointer to variance mixture fraction Z_v for each point of Z

    void ErrorMessage(const std::string error_message);
    int findingZResolution();

};

#endif


