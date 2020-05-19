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
	Enthalpy_Defect
Description
    Enthalpy_Defect is a step in the whole process, which reads the 1st
 	LUT to find deltaHSt. The class is a treatment for enthalpy defect 
	distribution
SourceFiles
    Enthalpy_Defect.C
\*---------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;


#ifndef Enthalpy_Defect_H
#define Enthalpy_Defect_H


class Enthalpy_Defect
{

public:

    Enthalpy_Defect();           // Constructor
	void Read();                        // Reading the LUT
    void Set_Path(const std::string ChiStName);
    
    void GetMeanValues(const double deltaH, double &extractedDeltaHSt);


private:
    std::string         _path_library;          // path to the LUT
    std::string         _name;                  // file name, in which deltaH distribution information is clarified 

    int     _ndeltaH;                                // number of mixture fraction deltaH
    std::vector<double>         _deltaH;             // deltaH values [J/kg]
    std::vector<double>         _deltaHSt;        // corresponding deltaHSt to the value of deltaH J/kg


    void ErrorMessage(const std::string error_message);
    int findingdeltaHResolution();

};

#endif


