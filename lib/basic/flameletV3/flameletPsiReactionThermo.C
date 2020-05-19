/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "flameletPsiReactionThermo.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(flameletPsiReactionThermo, 0);
    defineRunTimeSelectionTable(flameletPsiReactionThermo, fvMesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flameletPsiReactionThermo::flameletPsiReactionThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    psiThermo(mesh, phaseName),

    psi_
    (
        IOobject
        (
            phasePropertyName("thermo:psi"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, -2, 2, 0, 0)
    ),

    mu_
    (
        IOobject
        (
            phasePropertyName("thermo:mu"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    )


{
Info << " class flameletPsiReactionThermo is called\n" << endl;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::flameletPsiReactionThermo> Foam::flameletPsiReactionThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<flameletPsiReactionThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::flameletPsiReactionThermo::~flameletPsiReactionThermo()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::volScalarField& Foam::flameletPsiReactionThermo::Z()
{
    //notImplemented("flameletPsiReactionThermo::Z()");
    return const_cast<volScalarField&>(volScalarField::null());
}


const Foam::volScalarField& Foam::flameletPsiReactionThermo::Z() const
{
    //notImplemented("flameletPsiReactionThermo::Z() const");
    return volScalarField::null();
}

Foam::volScalarField& Foam::flameletPsiReactionThermo::Zvar()
{
    //notImplemented("flameletPsiReactionThermo::Zvar()");
    return const_cast<volScalarField&>(volScalarField::null());
}


const Foam::volScalarField& Foam::flameletPsiReactionThermo::Zvar() const
{
    //notImplemented("flameletPsiReactionThermo::Zvar() const");
    return volScalarField::null();
}


Foam::volScalarField& Foam::flameletPsiReactionThermo::chi_st()
{
   // notImplemented("flameletPsiReactionThermo::chi_st()");
    return const_cast<volScalarField&>(volScalarField::null());
}


const Foam::volScalarField& Foam::flameletPsiReactionThermo::chi_st() const
{
    //notImplemented("flameletPsiReactionThermo::chi_st() const");
    return volScalarField::null();
}


Foam::volScalarField& Foam::flameletPsiReactionThermo::H()
{
    //notImplemented("flameletPsiReactionThermo::H()");
    return const_cast<volScalarField&>(volScalarField::null());
}


const Foam::volScalarField& Foam::flameletPsiReactionThermo::H() const
{
    //notImplemented("flameletPsiReactionThermo::H() const");
    return volScalarField::null();
}

Foam::volScalarField& Foam::flameletPsiReactionThermo::as()
{
    Info << "\n as()  [1/s]Mean absorbtion coefficient of flameletPsiReactionThermo is called \n" << endl;
    //notImplemented("flameletPsiReactionThermo::as()");
    return const_cast<volScalarField&>(volScalarField::null());
}


const Foam::volScalarField& Foam::flameletPsiReactionThermo::as() const
{
    Info << "\n as() [1/s] Mean absorbtion coefficient of flameletPsiReactionThermo is called \n" << endl;
    //notImplemented("flameletPsiReactionThermo::as() const");
    return volScalarField::null();
}


Foam::tmp<Foam::volScalarField> Foam::flameletPsiReactionThermo::rho() const
{
    Info <<"\n rho [kg/m^3] in flameletPsiReactionThermo is updated by p_rgh_*psi_ \n" << endl;
    return p_*psi_;
}


Foam::tmp<Foam::scalarField> Foam::flameletPsiReactionThermo::rho(const label patchi) const
{
    Info <<"\n rho [kg/m^3] for patch in flameletPsiReactionThermo is updated p_rgh_.boundaryField()[patchi]*psi_.boundaryField()[patchi] \n" << endl;
    return p_.boundaryField()[patchi]*psi_.boundaryField()[patchi];
}


const Foam::volScalarField& Foam::flameletPsiReactionThermo::psi() const
{
    Info << "\n Compressibility psi_ [s^2/m^2] of flameletPsiReactionThermo is called \n" << endl;
    return psi_;
}


Foam::tmp<Foam::volScalarField> Foam::flameletPsiReactionThermo::mu() const
{
    Info << "\n Dynamic viscosity of mixture mu()[kg/m/s] of flameletPsiReactionThermo is called \n" << endl;
    return mu_;
}


Foam::tmp<Foam::scalarField> Foam::flameletPsiReactionThermo::mu(const label patchi) const
{
    Info << "\n Dynamic viscosity of mixture mu()[kg/m/s] for patch of flameletPsiReactionThermo is called \n" << endl;
    return mu_.boundaryField()[patchi];
}





// ************************************************************************* //
