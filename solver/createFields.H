
Info<< "Reading thermophysical properties\n" << endl;
autoPtr<flameletPsiReactionThermo> pThermo  // This function pThermo selecting thermo package in constant thermophysical
(                                // Reading all the information about PDF-Flamelet ... should look at in the log.application
    flameletPsiReactionThermo::New(mesh)    // This function first select thermalProperties in constant, run constructor in flameletPsiReactionThermo, thenpdfFlameletThermo
);

flameletPsiReactionThermo& thermo = pThermo(); 
thermo.validate(args.executable(), "h", "e");


// For new radiation model
basicMultiComponentMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.lookup("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species()
        << exit(FatalIOError);
}

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

// List of species in composition
forAll(Y,i)
{
Info << Y[i].name()<< endl;
}
fields.add(thermo.he());

#include "readAdditionalThermo.H"


Info<< "Creating field rho\n" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

Info<< "Preparing scalar fields from flameletPsiReactionThermo\n"
    << "  + thermo.p()\n"
    << "  + thermo.Z()\n"
    << "  + thermo.Zvar()\n"
    << "  + thermo.chi_st()\n"
    << "  + thermo.H()\n"
    << "  + thermo.as()\n"
    << endl;

volScalarField& p = thermo.p();
volScalarField& Z = thermo.Z();
volScalarField& Zvar = thermo.Zvar();
volScalarField& chi_st = thermo.chi_st();               // Not in "0" folder
volScalarField& H = thermo.H();
volScalarField& as = thermo.as();

Info<< "\nReading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField chi_stInstan
(
    IOobject
    (
        "chi_stInstan",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("chi_stInstan",dimensionSet(0,0,-1,0,0,0,0) , 0.0)
);

#include "compressibleCreatePhi.H"

#include "createMRF.H"


// Info<< "Creating turbulence~/OpenFOAM/leminh-dev/applications/line_burner_fireFoam-2/XO2-021/constant model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

const compressible::LESModel& sgs = turbulence->mesh().template lookupObject<compressible::LESModel>(turbulenceModel::propertiesName);



// Field dQ for radiation
volScalarField dQ
(
    IOobject
    (
        "dQ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dQ", dimEnergy/dimTime/dimVolume, 0.0) // kvm
);
#include "flameletControl.H"
#include "createRadiationModel.H"


Info<< "Preparing field Qrad (radiative heat transfer)\n" << endl;
volScalarField Qrad
(
    IOobject
    (
        "Qrad",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("qradiation", dimensionSet(1,-1,-3,0,0,0,0) , 0.0) // W/m^3
);

dimensionedScalar sigmaSB
    (
        "stefan_bolzmann",
        dimensionSet(1,0,-3,-4,0,0,0),
        5.6704e-8
    );






#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"
#include "readpRef.H"

volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

thermo.correct(); // kvm
rho = thermo.rho(); // kvm

mesh.setFluxRequired(p_rgh.name());

#include "phrghEqn.H"





Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

volScalarField dQrad
(
 IOobject
 (
  "dQrad",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("dQrad", dimEnergy/dimTime/dimVolume,0.0)

 );

volScalarField dQradabs
(
 IOobject
 (
  "dQradabs",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("dQradabs", dimEnergy/dimTime/dimVolume,0.0)

 );

volScalarField dQradem
(
 IOobject
 (
  "dQradem",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("dQradem", dimEnergy/dimTime/dimVolume,0.0)

 );

volScalarField IntermF
(
 IOobject
 (
  "IntermF",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("IntermF", dimless,0.0)
 );


volScalarField ZgradCells
(
 IOobject
 (
  "ZgradCells",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("ZgradCells", dimensionSet(0,-1,0,0,0,0,0), 0.0)
 );

volScalarField deltaCells
(
 IOobject
 (
  "deltaCells",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("deltaCells", dimensionSet(0,1,0,0,0,0,0), 0.0)
 );

// volScalarField hcF
// (
//     IOobject
//     (
//         "hcF",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("hcF", dimEnergy/dimMass, 0.0)
// );

// volScalarField heF
// (
//     IOobject
//     (
//         "heF",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("heF", dimEnergy/dimMass, 0.0)
// );


