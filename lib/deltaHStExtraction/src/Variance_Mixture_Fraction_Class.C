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

#include "Variance_Mixture_Fraction_Class.H"

Variance_Mixture_Fraction_Class::Variance_Mixture_Fraction_Class()
{
    _path_library           = "[unassigned]";
};  

void Variance_Mixture_Fraction_Class::Read()
{
    int i;
    string tag;
    _name = _path_library + "LookUpTable.out";
    //cout << _name << endl;
    _nZv = findingZvResolution();
    ifstream fInput;
    fInput.open(_name.c_str(), ios::in);
    fInput >> tag;
    if (tag != "Zv")
        ErrorMessage("Expected: Zv, Found: " + tag);
    fInput >> tag;
    if (tag != "index")
        ErrorMessage("Expected: index, Found: " + tag);
    _Zv.resize(_nZv);             // Resize a vector of _Zv
    _Zv_name.resize(_nZv);        // Resize a vector of _Zv_name

    for (i = 0; i < _nZv; i++)
    {
        fInput >> _Zv[i];
        fInput >> _Zv_name[i];
        //cout << _Zv[i] << "  " << _Zv_name[i]<< endl;         // print out to check
    }
        
    _ChiSt = new Stoich_Scalar_Dissipation_Rate_Class[_nZv];         // Dynamic allocation for each Z

    for (i = 0; i < _nZv; i++)
    {   
        _ChiSt[i].Set_Path(_path_library + _Zv_name[i] + "/");
        _ChiSt[i].Read();
    }
    

}


void Variance_Mixture_Fraction_Class::ErrorMessage(const string error_message)
{
    cout << "Variance_Mixture_Fraction_Class Error" << endl;
    cout << "File name: " << _name << endl;
    cout << "Error message: " << error_message << endl;
    exit(-1);
}

int Variance_Mixture_Fraction_Class::findingZvResolution()
{
    int ZvRes = 0;
    string line;
    ifstream fInput;
    fInput.open(_name.c_str(), ios::in);

    while (getline(fInput, line))
        ++ZvRes;
    return ZvRes - 1; // The first line is for label, so the number of Zv is equal the number of lines minus 1
}

void Variance_Mixture_Fraction_Class::Set_Path(const string ZName)
{
    _path_library = ZName ;
}

void Variance_Mixture_Fraction_Class::GetMeanValues( const double Zv, const double ChiSt, const double deltaH, double &extractedDeltaHSt)
{
    if (Zv >= _Zv[_nZv-1])         // greater or equal than max Zv in the 1st LUT
    {   
        _ChiSt[_nZv-1].GetMeanValues(ChiSt, deltaH, extractedDeltaHSt);
        //cout <<"Zv = " <<_Zv[_nZv-1] << "; folder name = " << _Zv_name[_nZv-1] << endl;
        return;
    }
    else if (Zv <= _Zv[0])            // smaller or equal than min ChiSt in the 1st LUT
    {
        _ChiSt[0].GetMeanValues(ChiSt, deltaH, extractedDeltaHSt);
        //cout <<"Zv = " <<_Zv[0] << "; folder name = " << _Zv_name[0] << endl;
        return;
    }
    else
    {
        int k = 0;
        double ratio;
        double extractedDeltaHSt_1;
        double extractedDeltaHSt_2;

        //finding k
        for( int j = 0; j < _nZv ; j++)
            if (Zv == _Zv[j])
            {
                _ChiSt[j].GetMeanValues(ChiSt, deltaH, extractedDeltaHSt);
                //cout <<"Zv = " <<_Zv[j] << "; folder name = " << _Zv_name[j] << endl;
                return;
            }
            else if(Zv < _Zv[j])
            {
                //cout << j << endl;
                k = j;
                break;
            }
        _ChiSt[k-1].GetMeanValues(ChiSt, deltaH, extractedDeltaHSt_1);
		_ChiSt[k].GetMeanValues(ChiSt, deltaH, extractedDeltaHSt_2);
        ratio = (Zv - _Zv[k-1])/(_Zv[k] - _Zv[k-1]);

        extractedDeltaHSt = extractedDeltaHSt_1 + ratio*(extractedDeltaHSt_2-extractedDeltaHSt_1);
        return;
    }
	

}




