    // for outputing flame height
    // OFstream outFlameHeight("outFlameHeight_"+runTime.timeName());
    // outFlameHeight << "# Time    " << "FlameHeight  "  << endl;

    // volScalarField flameHeight
    // (
    //     IOobject
    //     (
    //         "flameHeight",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ
    //     ),
    //     mesh,
    //     dimensionedScalar("zero",dimLength,0.0)
    // );
    // volScalarField flameHeight2
    // (
    //     IOobject
    //     (
    //         "flameHeight2",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ
    //     ),
    //     mesh,
    //     dimensionedScalar("zero",dimLength,0.0)
    // );

    // const pointField& cellCentres = mesh.cellCentres();


    // // estimation of HRR

    // OFstream outMassFlow("massFlow_"+runTime.timeName());
    // outMassFlow << "# Time    " << "massFlowO2 "  << endl;
    
    // OFstream outHRR("HRR"+runTime.timeName());
    // outHRR << "# Time    " << "HRR "  << endl;


    // OFstream outheF("heF"+runTime.timeName());
    // outHRR << "# Time    " << "heF "  << endl;

/**************************** Added by Rui to calcualte sensible enthalpy flux at the boundary *******************************/
    
    // create fileds for sensible enthalpy

    volScalarField heF
	(
    	IOobject
    	(
	        "heF",
	        runTime.timeName(),
	        mesh,
	        IOobject::READ_IF_PRESENT,
	        IOobject::AUTO_WRITE
    	),
    	mesh,
    	dimensionedScalar("heF", dimEnergy/dimMass, 0.0)
	);

    // create fileds for surface integration of sensible enthalpy flux

    surfaceScalarField phiHs
    (
        IOobject
        (
            "phiHs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi*fvc::interpolate(heF) 
    );




