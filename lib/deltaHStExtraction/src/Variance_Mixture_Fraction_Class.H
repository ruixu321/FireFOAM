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
	Variance_Mixture_Fraction_Class
Description
    Variance_Mixture_Fraction_Class is a step in the whole process, which reads the 1st
 	LUT to find deltaHSt. The class is a treatment for variance mixture fraction normalised 
	distribution
SourceFiles
    Variance_Mixture_Fraction_Class.C
\*---------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;


#ifndef Variance_Mixture_Fraction_Class_H
#define Variance_Mixture_Fraction_Class_H

#include "Stoich_Scalar_Dissipation_Rate_Class.H"

class Variance_Mixture_Fraction_Class
{

public:

    Variance_Mixture_Fraction_Class();           // Constructor
	void Read();                        // Reading the LUT
    void Set_Path(const std::string ZName);
    void GetMeanValues(const double Zv, const double ChiSt, const double deltaH, double &extractedDeltaHSt);


private:
    std::string         _path_library;          // path to the LUT
    std::string         _name;                  // file name, in which Zv distribution information is clarified 

    int     _nZv;                                // number of mixture fraction Zv
    std::vector<double>         _Zv;             // Zv values
    std::vector<string>         _Zv_name;        // corresponding directory name to the value of Zv
    Stoich_Scalar_Dissipation_Rate_Class         *_ChiSt;        // Pointer to each point of Z

    void ErrorMessage(const std::string error_message);
    int findingZvResolution();


};

#endif


