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

#include "Mixture_Fraction_Class.H"

Mixture_Fraction_Class::Mixture_Fraction_Class()
{
    _path_library           = "DeltaH/";


};

void Mixture_Fraction_Class::Read()
{
    int i;
    string tag;
    _name = _path_library + "LookUpTable.out";
    _nZ = findingZResolution();
    ifstream fInput;
    fInput.open(_name.c_str(), ios::in);
    fInput >> tag;
    if (tag != "Z")
        ErrorMessage("Expected: Z, Found: " + tag);
    fInput >> tag;
    if (tag != "index")
        ErrorMessage("Expected: index, Found: " + tag);
    _Z.resize(_nZ);             // Resize a vector of _Z
    _Z_name.resize(_nZ);        // Resize a vector of _Z_name

    for (i = 0; i < _nZ; i++)
    {
        fInput >> _Z[i];
        fInput >> _Z_name[i];
        //cout << _Z[i] << "  " << _Z_name[i]<< endl;         // print out to check
    }

    _Zv = new Variance_Mixture_Fraction_Class[_nZ];         // Dynamic allocation for each Z
    
    for (i = 0; i < _nZ; i++)
    {
        _Zv[i].Set_Path(_path_library + _Z_name[i] + "/");
        _Zv[i].Read();
    }

}


void Mixture_Fraction_Class::ErrorMessage(const string error_message)
{
    cout << "Mixture_Fraction_Class Error" << endl;
    cout << "File name: " << _name << endl;
    cout << "Error message: " << error_message << endl;
    exit(-1);
}

int Mixture_Fraction_Class::findingZResolution()
{
    int ZRes = 0;
    string line;
    ifstream fInput;
    fInput.open(_name.c_str(), ios::in);

    while (getline(fInput, line))
        ++ZRes;
    return ZRes - 1; // The first line is for label, so the number of Z is equal the number of lines minus 1
}


void Mixture_Fraction_Class::GetMeanValues(const double Z, const double Zv, const double ChiSt, const double deltaH, double &extractedDeltaHSt)
{

    if (Z >= _Z[_nZ-1])         // greater or equal than max Z in the 1st LUT
    {  
        _Zv[_nZ-1].GetMeanValues(Zv, ChiSt, deltaH, extractedDeltaHSt);
        //cout <<"Z = " <<_Z[_nZ-1] << "; folder name = " << _Z_name[_nZ-1] << endl;
        return;
    }
    else if (Z <= _Z[0])            // smaller or equal than min ChiSt in the 1st LUT
    {
        _Zv[0].GetMeanValues(Zv, ChiSt, deltaH, extractedDeltaHSt);
        //cout <<"Z = " <<_Z[0] << "; folder name = " << _Z_name[0] << endl;
        return;
    }
    else
    {
        int k = 0;
        double ratio;
        double extractedDeltaHSt_1;
        double extractedDeltaHSt_2;

        //finding k
        for( int j = 0; j < _nZ ; j++)
            if (Z == _Z[j])
            {
                _Zv[j].GetMeanValues(Zv, ChiSt, deltaH, extractedDeltaHSt);
                //cout <<"Z = " <<_Z[j] << "; folder name = " << _Z_name[j] << endl;
                return;
            }
            else if(Z < _Z[j])
            {
                //cout << j << endl;
                k = j;
                break;
            }
        _Zv[k-1].GetMeanValues(Zv, ChiSt, deltaH, extractedDeltaHSt_1);
        _Zv[k].GetMeanValues(Zv, ChiSt, deltaH, extractedDeltaHSt_2);
        ratio = (Z - _Z[k-1])/(_Z[k] - _Z[k-1]);

        extractedDeltaHSt = extractedDeltaHSt_1 + ratio*(extractedDeltaHSt_2-extractedDeltaHSt_1);
        return;
    }
   
}



