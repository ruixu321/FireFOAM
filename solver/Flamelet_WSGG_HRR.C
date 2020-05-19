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

Application
    flameletNewRadiationModel

Description
    Transient solver for fires and turbulent diffusion flames with reacting
    particle clouds, surface film and pyrolysis modelling.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "flameletPsiReactionThermo.H"
#include "turbulentFluidThermoModel.H"
#include "radiationModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "thermoPhysicsTypes.H"
#include "IOmanip.H"
#include "OFstream.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "printVersion.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "infoFieldsOutput.H" // fm
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"
	
    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;
        #include "rhoEqn.H"
	    #include "variablesControl.H"
            // --- PIMPLE loop
            while (pimple.loop())
            {

                #include "UEqn.H"
		        #include "ZEqn.H"
		        #include "HEqn.H"
             
	 	        Info << "After U Z and H equations" << endl;
	        // --- Pressure corrector loop
                while (pimple.correct())
                {
                    #include "pEqn.H"
                }

                if (pimple.turbCorr())
                {
                    turbulence->correct();
                }
            }

        rho = thermo.rho();

        #include "infoOutput.H"
	   	
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"; // kvm

        std::ofstream file;
        file.open ("time_total.txt", std::ofstream::out | std::ofstream::app);
        file << runTime.timeName() << " " << runTime.elapsedCpuTime() << std::endl << "\n";
        file << nl << endl;
        file.close();

        if (runTime.writeTime()) // kvm
        {
            Info<< " +"; // kvm
        }
        Info<< nl << endl; // kvm

    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
