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

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>


using namespace std;

#include "Enthalpy_Defect.H"

Enthalpy_Defect::Enthalpy_Defect()
{
    _path_library           = "[unassigned]";
};  

void Enthalpy_Defect::Read()
{
    int i;
    string tag;
    _name = _path_library + "PDF.csv";
    //cout << _name << endl;
    _ndeltaH = findingdeltaHResolution();
    //cout << "deltaH Resolution: " << _ndeltaH <<endl;
    ifstream fInput;
    fInput.open(_name.c_str(), ios::in);
	getline(fInput, tag, ',');
	if (tag != "DeltaH [J/kg]")
        ErrorMessage("Expected: DeltaH [J/kg], Found: " + tag);

    getline(fInput, tag);
    if (tag != "DeltaHst [J/kg]")
        ErrorMessage("Expected: DeltaHst [J/kg], Found: " + tag);

    _deltaH.resize(_ndeltaH);             // Resize a vector of _deltaH
    _deltaHSt.resize(_ndeltaH);        // Resize a vector of _deltaHSt

    for (i = 0; i < _ndeltaH; i++)
    {
		getline(fInput, tag, ',');
		_deltaH[i] = atof(tag.c_str());
		getline(fInput, tag);
		_deltaHSt[i] = atof(tag.c_str());
        //cout << _deltaH[i] << "  " << _deltaHSt[i]<< endl;         // print out to check
    }


}


void Enthalpy_Defect::ErrorMessage(const string error_message)
{
    cout << "Enthalpy_Defect Error" << endl;
    cout << "File name: " << _name << endl;
    cout << "Error message: " << error_message << endl;
    exit(-1);
}

int Enthalpy_Defect::findingdeltaHResolution()
{
    int deltaHRes = 0;
    string line;
    ifstream fInput;
    fInput.open(_name.c_str(), ios::in);

    while (getline(fInput, line))
        ++deltaHRes;
    return deltaHRes - 1; // The first line is for label, so the number of deltaH is equal the number of lines minus 1
}

void Enthalpy_Defect::Set_Path(const string ChiStName)
{
    _path_library = ChiStName ;
}


void Enthalpy_Defect::GetMeanValues( const double deltaH, double &extractedDeltaHSt)
{
    if ( _deltaH[0] == 0.0 )          // where Z = 0 or, Z =1, deltaH = 0, so deltaHSt should be 0. Rui 06182019 when deltaH ~ 0 , Hst = 0
    {
        extractedDeltaHSt = 0.0;
    }
	else if (deltaH >= _deltaH[_ndeltaH - 1]) // greater or equal than max deltaH in the 1st LUT
	{
		extractedDeltaHSt = _deltaHSt[_ndeltaH - 1];
		return;
	}
	else if (deltaH <= _deltaH[0])		// smaller or equal than min deltaH in the 1stLUT
	{
		extractedDeltaHSt = _deltaHSt[0];
		return;
	}
	else 
	{
		int k = 0;
		double ratio;
		
		//finding k
		for( int j = 0; j < _ndeltaH ; j++)
			if (deltaH == _deltaH[j])
			{
				extractedDeltaHSt	=	_deltaHSt[j];
				return;
			}
			else if(deltaH < _deltaH[j])
			{	
				//cout << j << endl;
				k = j;
				break;
			}
		ratio = (deltaH - _deltaH[k-1])/(_deltaH[k] - _deltaH[k-1]);
		extractedDeltaHSt = _deltaHSt[k-1] + ratio*(_deltaHSt[k] - _deltaHSt[k-1]);
	}	
	

}




