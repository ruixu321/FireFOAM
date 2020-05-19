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

#include "Stoich_Scalar_Dissipation_Rate_Class.H"

Stoich_Scalar_Dissipation_Rate_Class::Stoich_Scalar_Dissipation_Rate_Class()
{
    _path_library           = "[unassigned]";
};  

void Stoich_Scalar_Dissipation_Rate_Class::Read()
{
    int i;
    string tag;
    _name = _path_library + "LookUpTable.out";
    //cout << _name << endl;
    _nChiSt = findingChiStResolution();
    //cout << "ChiSt Resolution: " << _nChiSt <<endl;
    ifstream fInput;
    fInput.open(_name.c_str(), ios::in);
    fInput >> tag;
    if (tag != "chist")
        ErrorMessage("Expected: ChiSt, Found: " + tag);
    fInput >> tag;
    if (tag != "index")
        ErrorMessage("Expected: index, Found: " + tag);
    _ChiSt.resize(_nChiSt);             // Resize a vector of _ChiSt
    _ChiSt_name.resize(_nChiSt);        // Resize a vector of _ChiSt_name

    for (i = 0; i < _nChiSt; i++)
    {
        fInput >> _ChiSt[i];
        fInput >> _ChiSt_name[i];
        //cout << _ChiSt[i] << "  " << _ChiSt_name[i]<< endl;         // print out to check
    }

    _deltaH = new Enthalpy_Defect[_nChiSt];         // Dynamic allocation for each Z
    
    for (i = 0; i < _nChiSt; i++)
    {
        _deltaH[i].Set_Path(_path_library + _ChiSt_name[i] + "/");
        _deltaH[i].Read();
    }


}


void Stoich_Scalar_Dissipation_Rate_Class::ErrorMessage(const string error_message)
{
    cout << "Stoich_Scalar_Dissipation_Rate_Class Error" << endl;
    cout << "File name: " << _name << endl;
    cout << "Error message: " << error_message << endl;
    exit(-1);
}

int Stoich_Scalar_Dissipation_Rate_Class::findingChiStResolution()
{
    int ChiStRes = 0;
    string line;
    ifstream fInput;
    fInput.open(_name.c_str(), ios::in);

    while (getline(fInput, line))
        ++ChiStRes;
    return ChiStRes - 1; // The first line is for label, so the number of ChiSt is equal the number of lines minus 1
}

void Stoich_Scalar_Dissipation_Rate_Class::Set_Path(const string ZvName)
{
    _path_library = ZvName ;
}

void Stoich_Scalar_Dissipation_Rate_Class::GetMeanValues( const double ChiSt, const double deltaH, double &extractedDeltaHSt)
{
    if (ChiSt >= _ChiSt[_nChiSt-1])         // greater or equal than max ChiSt in the 1st LUT
    {    
        _deltaH[_nChiSt-1].GetMeanValues(deltaH, extractedDeltaHSt);
        //cout <<"ChiSt = " <<_ChiSt[_nChiSt-1] << "; folder name = " << _ChiSt_name[_nChiSt-1] << endl;
        return;
    }
    else if (ChiSt <= _ChiSt[0])            // smaller or equal than min ChiSt in the 1st LUT
    {
        _deltaH[0].GetMeanValues(deltaH, extractedDeltaHSt);
        //cout <<"ChiSt = " <<_ChiSt[0] << "; folder name = " << _ChiSt_name[0] << endl;
        return;
    }
    else
    {
        int k = 0;
        double ratio;
        double extractedDeltaHSt_1;
        double extractedDeltaHSt_2;
        
        //finding k
        for( int j = 0; j < _nChiSt ; j++)
            if (ChiSt == _ChiSt[j])
            {
                _deltaH[j].GetMeanValues(deltaH, extractedDeltaHSt);
				//cout <<"ChiSt = " <<_ChiSt[j] << "; folder name = " << _ChiSt_name[j] << endl;
                return;
            }
            else if(ChiSt < _ChiSt[j])
            {
                //cout << j << endl;
                k = j;
                break;
            }
		_deltaH[k-1].GetMeanValues(deltaH, extractedDeltaHSt_1);
		_deltaH[k].GetMeanValues(deltaH, extractedDeltaHSt_2);
		ratio = (ChiSt - _ChiSt[k-1])/(_ChiSt[k] - _ChiSt[k-1]);
		
		extractedDeltaHSt = extractedDeltaHSt_1 + ratio*(extractedDeltaHSt_2-extractedDeltaHSt_1);
		return;	
    }
}





