
// Rui 07/23/2019
double time_rad = runTime.elapsedCpuTime();

// Info<< "ExecutionTime_Start_Radiation= " << runTime.elapsedCpuTime() << " s";
radiation->correct(); 
// Info<< "ExecutionTime_End_Radiation= " << runTime.elapsedCpuTime() << " s";

time_rad = runTime.elapsedCpuTime() - time_rad;

std::ofstream file;
file.open ("time_rad.txt", std::ofstream::out | std::ofstream::app);
file << runTime.timeName() << " " << time_rad << std::endl << "\n";
file << nl << endl;
file.close();


// Minh's codes to use absorption coefficient and emission with TRI treatment from flamelet LUT
const volScalarField& SumMHFl = mesh.lookupObject<volScalarField>("SumMHFl");
// const volScalarField& ProdProgVarFl = mesh.lookupObject<volScalarField>("ProdProgVarFl");


const volScalarField& EmTRI1 = mesh.lookupObject<volScalarField>("Em1TRI");
const volScalarField& EmTRI2 = mesh.lookupObject<volScalarField>("Em2TRI");
const volScalarField& EmTRI3 = mesh.lookupObject<volScalarField>("Em3TRI");
const volScalarField& EmTRI4 = mesh.lookupObject<volScalarField>("Em4TRI");

const volScalarField EmTRIMul4Pi = EmTRI1*4.0*constant::mathematical::pi + EmTRI2*4.0*constant::mathematical::pi + 
                                    EmTRI3*4.0*constant::mathematical::pi + EmTRI4*4.0*constant::mathematical::pi ;

const volScalarField::Internal rad = radiation->Ru() - EmTRIMul4Pi;  // Absorption - emission (total)
volScalarField::Internal& radField = dQrad.ref();
radField = rad;

const volScalarField::Internal radRu = radiation->Ru();    // Absorption
volScalarField::Internal& radFieldRu = dQradabs.ref();
radFieldRu = radRu;

const volScalarField::Internal radRp = EmTRIMul4Pi;       // emission
volScalarField::Internal& radFieldRp = dQradem.ref();
radFieldRp = radRp;

{
    volScalarField dEff
    (
        "dEff",
        turbulence->alphaEff()
    );

    if (!constD)
    {
        dEff -= turbulence->alpha()*(1-1.0/lewisNo);
    }
    else
    {
        dEff = dEff - turbulence->alpha() + DM*rho;
    }
    
    if (adiabaticMode == true)
        {
        Info << "In the adiabatic mode" << endl;
            if (radiationMode == true)
            {
                    Qrad = 4. * sigmaSB * as * (pow(thermo.T(),4) - pow(Tenv,4));

            }

            else if(radiationMode == false && switchQoff == false)
            {
                    Qrad = 0 * 4. * sigmaSB * as * (pow(thermo.T(),4) - pow(Tenv,4));
                    switchQoff = true;
            }


            fvScalarMatrix hEqn
            (
            fvm::ddt(rho, H)	
              + fvm::div(phi, H)
              - fvm::laplacian(dEff, H)      
            );
            hEqn.relax();
            hEqn.solve(mesh.solver("H"));
        }

    else
        {
            Info << "In the non-adiabatic mode" << endl;
            fvScalarMatrix hEqn
            (
            fvm::ddt(rho, H)	
              + fvm::div(phi, H)
              - fvm::laplacian(dEff, H)  - dQrad
            );

            hEqn.relax();
            hEqn.solve(mesh.solver("H"));
        }
    thermo.correct();
    forAll(Y, i)
       {
       Info << setw(8) << Y[i].name() << "\tmin/ave/max "<< " =\t"
            << setw(12) << setprecision(5) << gMin(Y[i]) << "\t"
            << setw(12) << setprecision(5) << gAverage(Y[i]) << "\t"
            << setw(12) << setprecision(5) << gMax(Y[i]) << nl;
       }
}
