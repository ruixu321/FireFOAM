/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::pdfFlameletPsiReactionThermo

Description
    Energy for a mixture based on density

SourceFiles
    pdfFlameletPsiReactionThermo.C

\*---------------------------------------------------------------------------*/

#ifndef pdfFlameletPsiReactionThermo_H
#define pdfFlameletPsiReactionThermo_H

#include "flameletPsiReactionThermo.H"
#include "heThermo.H"	//Enthalpy/Internal energy for a mixture
#include "OpenSMOKE_PDF_NonAdiabaticFlamelet_Library.H"
#include "Mixture_Fraction_Class.H"         // Class for deltaHSt extraction

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pdfFlameletPsiReactionThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicPsiThermo, class MixtureType>
class pdfFlameletPsiReactionThermo
:
    public heThermo<BasicPsiThermo, MixtureType>  // We inherit heThermo class with specialization template BasicPsiThermo and MixtureType
{
    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();

        //- Update hte thermo variables by extracting from LUT
        void update();

        //- Update the species mass fraction by extracting from LUT
        void updateMassFractions();

        //- Error message
        void errorMessage(const string message);

        //- Info message
        void infoMessage() const;

        //- Construct as copy (not implemented)
        pdfFlameletPsiReactionThermo
        (
            const pdfFlameletPsiReactionThermo<BasicPsiThermo,
            MixtureType>&
        );


public:

    //- Runtime type information
    TypeName("pdfFlameletPsiReactionThermo");


    // Constructors

        //- Construct from mesh and phase name
        pdfFlameletPsiReactionThermo
        (
            const fvMesh&,
            const word& phaseName
        );


    //- Destructor
    virtual ~pdfFlameletPsiReactionThermo();


    // Member functions

        //- Update properties
        virtual void correct();

private:

    // Variables

        //- Mixture fraction
        volScalarField Z_;

        //- Variance of mixture fraction
        volScalarField Zvar_;

        //- Scalar dissipation rate
        volScalarField chi_st_;

        //- Enthalpy of mixture 
        volScalarField H_;

        //- Enthalpy defect due to non-adiabatic system
        volScalarField defect_;

        //- Enthalpy defect at stoichiometry due to non-adiabatic system
        volScalarField defectSt_;

        //- PLank absorbtions coeff (for radiation)
        volScalarField as_;

        //- Reynolds density
        volScalarField density_reynolds_;

        //- Dynamic viscosity (favre)
        volScalarField mu_favre_;

        //- Thermal diffusivity
        volScalarField alpha_favre_;
        
        volScalarField MW_;
        volScalarField Cp_;
        volScalarField lambda_;
       
        //- PLank absorbtions coeff (for radiation) (1/s) (reynolds)
        volScalarField kappa1_reynolds_;             // Added by Rui to have 5 bands WSGG
        volScalarField kappa2_reynolds_;
        volScalarField kappa3_reynolds_;
        volScalarField kappa4_reynolds_;
         
        //- Radiative directional emission (W/m3/sr)       (reynolds)
        volScalarField Em1_reynolds_;
        volScalarField Em2_reynolds_;
        volScalarField Em3_reynolds_;
        volScalarField Em4_reynolds_;

        //- Heat release rate (W/m3) (reynolds)
        volScalarField SumMH_reynolds_;

        //- Molar fraction of N2                 // Added by Rui
        volScalarField MolarFracN2_favre_;

        //- Molar fraction of O2
        volScalarField MolarFracO2_favre_;

        //- Production rate of Progress Variable       (reynolds)
        // volScalarField ProdProgVar_reynolds_;


        //- Species mass fraction
        PtrList<volScalarField> omega_;				// A templated 1D list of pointers to objects of type <T>

        //- Adiabatic enthalpy of fuel
        scalar HFuel;

        //- Adiabatic enthalpy of oxidizer
        scalar HOxidizer;

        //- Variable for set thermo variable extraction
        label counter;						// Have to know what is label, maybe we can print ...

        //- variable for set mass fraction extraction
        label counter_mass_fractions;

        //- User defined numbers for update
        label propertyUpdate;
        label massFractionsUpdate;

        //label nDirac;
        Switch adiabaticMode;
        Switch showFlamelet;
        Switch showFlameletLibrary;

        //- stuff (modifyed by Holzmann)
        List<int> patch_type_T;
        List<int> patch_type_H;
        List<int> patch_type_Z;


private:

    //- Flamelet class
    OpenSMOKE_PDF_NonAdiabaticFlamelet_Library flamelets_library;

    //- Class to extract defectSt (deltaHSt)
    Mixture_Fraction_Class deltaHStExtractionObject;

public:

    //- Mixture fraction

        //  Non-const access allowed for transport equations
        virtual volScalarField& Z()
        {
            Info << "\n Non-const access of Z() of pdfFlameletPsiReactionThermo is called \n" << endl;   
            return Z_;
        }

        //- Mixture fraction
        virtual const volScalarField& Z() const
        {
            Info << "\n const of Z() of pdfFlameletPsiReactionThermo is called \n" << endl;   
            return Z_;
        }

        //- Variance of Mixture fraction
        //  Non-const access allowed for transport equations
        virtual volScalarField& Zvar()
        {
            Info << "\n Non-const access of Zvar() of pdfFlameletPsiReactionThermo is called \n" << endl;   
            return Zvar_;
        }

        //- Variance of Mixture fraction
        virtual const volScalarField& Zvar() const
        {
            Info << "\n const access of Zvar() of pdfFlameletPsiReactionThermo is called \n" << endl;   
            return Zvar_;
        }

        //- Scalar dissipation rate [1/s]
        //  Non-const access allowed for transport equations
        virtual volScalarField& chi_st()
        {
            Info << "\n Non-const access of chi_st() of pdfFlameletPsiReactionThermo is called \n" << endl;   
            return chi_st_;
        }

        //- Scalar dissipation rate [1/s]
        virtual const volScalarField& chi_st() const
        {
            Info << "\n Const access of chi_st() of pdfFlameletPsiReactionThermo is called \n" << endl;   
            return chi_st_;
        }

        //- Enthalpy [J/kg]
        //  Non-const access allowed for transport equations
        virtual volScalarField& H()
        {
            Info << "\n Non-const access of H() of pdfFlameletPsiReactionThermo is called \n" << endl;   
            return H_;
        }

        //- Enthalpy [J/kg]
        virtual const volScalarField& H() const
        {
            Info << "\n Const access H() of pdfFlameletPsiReactionThermo is called \n" << endl;   
            return H_;
        }

        //- Absorption coefficient [1/m]
        //  Non-const access allowed for transport equations
        virtual volScalarField& as()
        {
	    Info << "\n as() Mean absorption coefficient of pdfFlameletPsiReactionThermo is called \n" << endl;	
            return as_;
        }

        //- Absorption coefficient [1/m]
        virtual const volScalarField& as() const
        {
            Info << "\n const as()  Mean absorption coefficient of pdfFlameletPsiReactionThermo is called \n" << endl;   
            return as_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "pdfFlameletPsiReactionThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
