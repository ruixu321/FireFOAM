/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::flameletPsiReactionThermo

Description
    Foam::flameletPsiReactionThermo

SourceFiles
    flameletPsiReactionThermo.C

\*---------------------------------------------------------------------------*/

#ifndef flameletPsiReactionThermo_H
#define flameletPsiReactionThermo_H

#include "psiThermo.H"
#include "basicSpecieMixture.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class flameletPsiReactionThermo Declaration
\*---------------------------------------------------------------------------*/

class flameletPsiReactionThermo
:
    public psiThermo
{

protected:

    // Protected data

        //- Compressibility [s^2/m^2]
        volScalarField psi_;

        //- Dynamic viscosity [kg/m/s]
        volScalarField mu_;



public:

    //- Runtime type information
    TypeName("flameletPsiReactionThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        flameletPsiReactionThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Constructors

        //- Construct from mesh and phase name
        flameletPsiReactionThermo
        (
            const fvMesh&,
            const word& phaseName
        );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<flameletPsiReactionThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~flameletPsiReactionThermo();


    // Member functions

        //- Return the composition of the multi-component mixture
        virtual basicSpecieMixture& composition() = 0;

        //- Return the composition of the multi-component mixture
        virtual const basicSpecieMixture& composition() const = 0;



    // Member functions Minh

        //- Fields comming with pdfThermo


            //- Mixture fraction [-] non-const
            virtual volScalarField& Z();

            //- Mixture fraction [-] const
            virtual const volScalarField& Z() const;


            //- Varianz of mixture fraction [-] non-const
            virtual volScalarField& Zvar();

            //- Varianz of mixture fraction [-] const
            virtual const volScalarField& Zvar() const;

            //- Stochiometric mixture fraction [1/s] non-const
            virtual volScalarField& chi_st();

            //- Stochiometric mixture fraction [1/s] const
            virtual const volScalarField& chi_st() const;

            //- Enhalpy [J/kg] non-const
            virtual volScalarField& H();

            //- Enhalpy [J/kg] const
            virtual const volScalarField& H() const;

            //- Mean absorbtion coefficient [1/m] non-const
            virtual volScalarField& as();

            //- Mean absorbtion coefficient [1/m] const
            virtual const volScalarField& as() const;


        // Fields derived from thermodynamic state variables

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const;

            //- Compressibility [s^2/m^2]
            virtual const volScalarField& psi() const;


        // Access to transport state variables

            //- Dynamic viscosity of mixture [kg/m/s]
            virtual tmp<volScalarField> mu() const;

            //- Dynamic viscosity of mixture for patch [kg/m/s]
            virtual tmp<scalarField> mu(const label patchi) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
