// Estimation of flame height


// const scalar ZSt = 0.055; // For the methane fuel

// if (mag(g).value() > 0.0) 
// {    
//     flameHeight = dimensionedScalar("zero",dimLength,0.0);
//     flameHeight2 = dimensionedScalar("zero",dimLength,0.0);

//     forAll (Z.internalField(),cellI)
//     {
//         if (Z[cellI] >= ZSt)
//         { 
//             flameHeight[cellI] = mag( cellCentres[cellI] & (g/mag(g)).value() );
//         }
//         if (Z[cellI] >= 0.2*ZSt)
//         { 
//             flameHeight2[cellI] = mag( cellCentres[cellI] & (g/mag(g)).value() );
//         }
//     }

//     outFlameHeight <<  runTime.value() << "  "     
//             << max(flameHeight).value() << "  "
//             << max(flameHeight2).value() << "  "
//             << endl;
// }


// estimation of HRR

// const volScalarField& O2 = mesh.lookupObject<volScalarField>("O2");
// const volScalarField::Boundary& O2Bf = O2.boundaryField();

// scalar TotalsumFlux = 0.0;
// scalar TotalsumFluxcoflow = 0.0;

// forAll(O2Bf, patchi){
// 	word BCtype = O2Bf.types()[patchi];

// 	if(BCtype != "processor"){
// 		word patchlabelName = mesh.boundaryMesh().names()[patchi];

// 		const scalarField & phiw = phi.boundaryField()[patchi];
// 		const scalarField & O2w = O2.boundaryField()[patchi];

// 		// Initialise the sumFlux variable

// 		scalar sumFlux = 0.0;

// 		// O2 mass flow rate

// 		const scalarField wO2 = phiw*O2w;

// 		// Sum all of the face values and assign it to sumFlux
// 		sumFlux = Foam::sum(wO2);


// 		//Info << "Mass flow rate on patch " << patchlabelName << " : " << sumFlux <<  " kg/s" << endl;

// 		TotalsumFlux += sumFlux;
		
// 		// mass flux of O2 at coflow for checking
// 		if (patchlabelName == "coflow")
// 		TotalsumFluxcoflow =  sumFlux;

// 	}
// }

// reduce(TotalsumFlux, sumOp<scalar>() );
// reduce(TotalsumFluxcoflow, sumOp<scalar>() );

// Info << "Mass flow rate O2 at coflow : " << TotalsumFluxcoflow <<  " kg/s " << endl;
// Info << "Mass flow rate O2 : " << TotalsumFlux <<  " kg/s " << endl;

// scalar HRR = 12.5E6 * TotalsumFlux;



// outMassFlow <<  runTime.value() << "  "     
// 		<< TotalsumFlux 
// 		<< endl;

// outHRR << runTime.value() << " " 
//        << HRR  
//        << endl;

//Estimation of HRR by he
// heF = H - thermo.hc();

// const volScalarField::Boundary& heFBf = heF.boundaryField();

// scalar TotalheF = 0.0;
// scalar TotalhcF = 0.0;

// forAll(heFBf, patchi){
// 	word BCtype = heFBf.types()[patchi];

// 	if(BCtype != "processor"){
// 		word patchlabelName = mesh.boundaryMesh().names()[patchi];

// 		const scalarField & phiw = phi.boundaryField()[patchi];
// 		const scalarField & heFw = heF.boundaryField()[patchi];
// 		const volScalarField hcF = thermo.hc();
// 		const scalarField & hcFw = hcF.boundaryField()[patchi];

// 		// Initialise the sumFlux variable

// 		scalar sumFlux = 0.0;
// 		scalar sumFluxhc = 0.0;

// 		const scalarField wheF = phiw*heFw;
// 		const scalarField whcF = phiw*hcFw;

// 		// Sum all of the face values and assign it to sumFlux
// 		sumFlux = Foam::sum(wheF);
// 		sumFluxhc = Foam::sum(whcF);


// 		//Info << "Mass flow rate on patch " << patchlabelName << " : " << sumFlux <<  " kg/s" << endl;

// 		TotalheF += sumFlux;
// 		TotalhcF += sumFluxhc;

// 	}
// }

// reduce(TotalheF, sumOp<scalar>() );
// reduce(TotalhcF, sumOp<scalar>() );


// outheF << runTime.value() << " " 
//        << TotalheF << " "
//        << TotalhcF  
//        << endl;




/**************************** Added by Rui to calcualte sensible enthalpy flux at the boundary *******************************/

tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

heF = H - thermo.hc();

phiHs = phi*mvConvection->interpolate(phi,heF)
         - fvc::interpolate(turbulence->alphaEff())*fvc::snGrad(heF)*mesh.magSf();