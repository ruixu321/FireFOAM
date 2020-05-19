Info<< "Reading flameletProperties dictionary\n" << endl;       // flamelet
IOdictionary flameletsProperties
(
    IOobject
    (
        "flameletsProperties",
        U.time().constant(),
        U.db(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
dimensionedScalar initialMass = fvc::domainIntegrate(rho);
//- Monitoring
Switch monitoring(flameletsProperties.lookup("monitoring"));

//- Semi-empirical constants in mixture fraction variance equation
scalar Cg = readScalar(flameletsProperties.lookup("Cg"));
scalar Cd = readScalar(flameletsProperties.lookup("Cd"));


//- Semi-empirical constant for scalar dissipation rate evaluation
scalar C_chi = readScalar(flameletsProperties.lookup("C_chi"));

//- Semi-empirical constant for variance of mixture fraction
scalar C_Z = readScalar(flameletsProperties.lookup("C_Z"));

//- Prandtl number for mixture fraction, variance of mixture fraction
//  and enthalpy diffusion terms
dimensionedScalar sigmat
    (
        "sigmat",
        dimensionSet(0,0,0,0,0,0,0),
        readScalar(flameletsProperties.lookup("sigmat"))
    );
//- Minimum values (to be used for scalar dissipation rate evaluation)
dimensionedScalar k_small(flameletsProperties.lookup("k_small"));
dimensionedScalar epsilon_small
    (
        flameletsProperties.lookup("epsilon_small")
    );

//- Environment temperature (non adiabatic consitions)
dimensionedScalar Tenv (flameletsProperties.lookup("Tenv"));

//- Variance of mixture fraction:
//      on = transport equation 
//      off = algebraic equation
Switch varianceEquation(flameletsProperties.lookup("varianceEquation"));

//- Adiabatic simulation:
//      on = adiabatic
//      off = non adiabatic
Switch adiabaticMode(flameletsProperties.lookup("adiabaticMode"));
Switch radiationMode(flameletsProperties.lookup("radiationMode"));
bool switchQoff = false;

