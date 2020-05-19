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

\*---------------------------------------------------------------------------*/

#include "pdfFlameletPsiReactionThermo.H"
#include "Time.H"
#include "fixedValueFvPatchFields.H" // This boundary condition supplies a fixed value constraint, and is the base class for a number of other boundary conditions

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::pdfFlameletPsiReactionThermo<BasicPsiThermo, MixtureType>::calculate()
{
    const scalarField& pCells = this->p_.internalField();

    const scalarField& HCells = this->H_.internalField();

    const scalarField& ZCells = this->Z_.internalField();

    const scalarField& RhoReynolds = this->density_reynolds_.internalField();

    const scalarField& muFavre = this->mu_favre_.internalField();

    const scalarField& alphaFavre = this->alpha_favre_.internalField();

    const scalarField& TCells = this->T_.internalField(); //Minh

    scalarField& psiCells = this->psi_.ref();

    scalarField& muCells = this->mu_.ref();

    scalarField& alphaCells = this->alpha_.ref();

    scalarField& defectCells = this->defect_.ref();


    forAll(ZCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        //psiCells[celli] = RhoReynolds[celli]/pCells[celli];
        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);   //Minh

        //muCells[celli] = muFavre[celli];
	    muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);

        //alphaCells[celli] = alphaFavre[celli];
	    alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);

        defectCells[celli] = 
            HCells[celli] - (HOxidizer+ZCells[celli]*(HFuel-HOxidizer));
        
    }

    // Boundaries
    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];

        const fvPatchScalarField& pH = this->H_.boundaryField()[patchi];

        const fvPatchScalarField& pZ = this->Z_.boundaryField()[patchi];

        const fvPatchScalarField& pRhoReynolds
            = this->density_reynolds_.boundaryField()[patchi];

        const fvPatchScalarField& pmuFavre
            = this->mu_favre_.boundaryField()[patchi];

        const fvPatchScalarField& palphaFavre
            = this->alpha_favre_.boundaryField()[patchi];

        fvPatchScalarField& pT = this->T_.boundaryFieldRef()[patchi];

        fvPatchScalarField& ppsi = this->psi_.boundaryFieldRef()[patchi];

        fvPatchScalarField& pdefect = this->defect_.boundaryFieldRef()[patchi];

        fvPatchScalarField& pmu = this->mu_.boundaryFieldRef()[patchi];

        fvPatchScalarField& palpha = this->alpha_.boundaryFieldRef()[patchi];


        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);		//Minh

                //ppsi[facei] = pRhoReynolds[facei]/pp[facei];
                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]); //Minh

                //pmu[facei] = pmuFavre[facei];
		        pmu[facei] = mixture_.mu(pp[facei], pT[facei]); //Minh


                //palpha[facei] = palphaFavre[facei];
		        palpha[facei] = mixture_.alphah(pp[facei], pT[facei]); //Minh

                pdefect[facei]
                    = pH[facei] - (HOxidizer+pZ[facei]*(HFuel-HOxidizer));
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);			//Minh

                //ppsi[facei] = pRhoReynolds[facei]/pp[facei];
                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);		//Minh

                //pmu[facei] = pmuFavre[facei];
		        pmu[facei] = mixture_.mu(pp[facei], pT[facei]);               //Minh
		
                //palpha[facei] = palphaFavre[facei];
		        palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);               //Minh

                pdefect[facei]
                    = pH[facei] - (HOxidizer+pZ[facei]*(HFuel-HOxidizer));
            }
        }
    }
}


template<class BasicFlameletThermo, class MixtureType>
void Foam::pdfFlameletPsiReactionThermo<BasicFlameletThermo, MixtureType>::update()
{
    std::vector<scalar> extracted(20); //Old:7

    scalar Zvar_normalized = 0.;

    scalar defect = 0.;

    scalar defectSt = 0.0;      //Minh

    const scalarField& Z = this->Z_.internalField();

    const scalarField& Zvar = this->Zvar_.internalField();

    const scalarField& chi_st = this->chi_st_.internalField();

    const scalarField& HCells = this->H_.internalField();

    scalarField& defectStCells = this->defectSt_.ref();    

    scalarField& TCells = this->T_.ref();

    scalarField& RhoCells = this->density_reynolds_.ref();

    scalarField& asCells = this->as_.ref();

    scalarField& MWCells = this->MW_.ref(); // Rui
    scalarField& CpCells = this->Cp_.ref(); // Rui
    scalarField& lambdaCells = this->lambda_.ref(); // Rui
    
    
    scalarField& kappa1Cells = this->kappa1_reynolds_.ref();      //Rui
    scalarField& kappa2Cells = this->kappa2_reynolds_.ref();
    scalarField& kappa3Cells = this->kappa3_reynolds_.ref();
    scalarField& kappa4Cells = this->kappa4_reynolds_.ref();

    scalarField& Em1Cells = this->Em1_reynolds_.ref();      //Rui
    scalarField& Em2Cells = this->Em2_reynolds_.ref();
    scalarField& Em3Cells = this->Em3_reynolds_.ref();
    scalarField& Em4Cells = this->Em4_reynolds_.ref();

    scalarField& SumMHCells = this->SumMH_reynolds_.ref();      //Minh

    // scalarField& ProdProgVarCells = this->ProdProgVar_reynolds_.ref();

    scalarField& MolarFracN2Cells = this->MolarFracN2_favre_.ref();    // Rui
    scalarField& MolarFracO2Cells = this->MolarFracO2_favre_.ref();   // Rui


    scalarField& muCells = this->mu_favre_.ref();
    
    scalarField& alphaCells = this->alpha_favre_.ref();

    scalar small_eps = 1.e-6;

    scalar small_chi_st = 1.e-8;

    //- Internal cells
    forAll(Z, celli)
    {
        double max_chi = max(small_chi_st,chi_st[celli]);

        if (adiabaticMode == false)
        {
            defect = HCells[celli] - (HOxidizer+Z[celli]*(HFuel-HOxidizer));
            
        }

        //- Pure oxidizer
        if (Z[celli]<=small_eps)
        {
            deltaHStExtractionObject.GetMeanValues(0.,0.,max_chi,defect, defectSt); //Minh
            flamelets_library.GetMeanValues
            (
                0.,
                0.,
                max_chi,
                defectSt,                   // Old: defect, Minh: defectSt
                extracted
            );
        }

        //- Pure fuel
        else if (Z[celli]>=(1.-small_eps))
        {
            deltaHStExtractionObject.GetMeanValues(1.,0.,max_chi,defect, defectSt); //Minh

            flamelets_library.GetMeanValues
            (
                1.,
                0.,
                max_chi,
                defectSt,                   // Old: defect, Minh: defectSt
                extracted
            );
        }

        //- Mixture
        else
        {
            Zvar_normalized = Zvar[celli] / (Z[celli]*(1.-Z[celli]));

            if (Zvar_normalized >= 0.98)
            {
                deltaHStExtractionObject.GetMeanValues(Z[celli], 0.98,max_chi,defect, defectSt); //Minh        

                flamelets_library.GetMeanValues
                (
                    Z[celli],
                    0.98,
                    max_chi,
                    defectSt,           // Old: defect, Minh: defectSt
                    extracted
                );
            }
            else if (Zvar_normalized < 0.)
            {
                deltaHStExtractionObject.GetMeanValues(Z[celli], 0.0,max_chi,defect, defectSt); //Minh    

                flamelets_library.GetMeanValues
                (
                    Z[celli],
                    0.00,
                    max_chi,
                    defectSt,               // Old: defect, Minh: defectSt
                    extracted
                );
            }
            else
            {
                deltaHStExtractionObject.GetMeanValues(Z[celli], Zvar_normalized,max_chi,defect, defectSt); //Minh

                flamelets_library.GetMeanValues
                (
                    Z[celli],
                    Zvar_normalized,
                    max_chi,
                    defectSt,               // Old: defect, Minh: defectSt
                    extracted
                );
            }
        }
    
        defectStCells[celli] = defectSt;

        TCells[celli] = extracted[1];

        RhoCells[celli] = extracted[2];

        asCells[celli] = extracted[3];

        MWCells[celli] = extracted[4]  / 1000.0;  // g/mole to kg/mole

        CpCells[celli] = extracted[5];

        lambdaCells[celli] = extracted[6];

        muCells[celli] = extracted[7];

        alphaCells[celli] = extracted[8];

        kappa1Cells[celli] = extracted[9];    //Minh
        kappa2Cells[celli] = extracted[10];
        kappa3Cells[celli] = extracted[11];
        kappa4Cells[celli] = extracted[12];

        Em1Cells[celli] = extracted[13];      //Minh
        Em2Cells[celli] = extracted[14];
        Em3Cells[celli] = extracted[15];
        Em4Cells[celli] = extracted[16];

        SumMHCells[celli] = extracted[17];    //Minh

        MolarFracN2Cells[celli] = extracted[18];    //Rui
        MolarFracO2Cells[celli] = extracted[19];    //Rui

        // ProdProgVarCells[celli] = extracted[16];      //Minh

	
	// Relaxation until 5 second of properties T, rho, as, mu, alpha
	scalar x = this->time().value();
	//TCells[celli]	 = TCells[celli]*(0.5+0.5*Foam::tanh((x-0.500)/0.250)) + 294.75*(1-(0.5+0.5*Foam::tanh((x-0.500)/0.250)));
	//RhoCells[celli]	 = RhoCells[celli]*(0.5+0.5*Foam::tanh((x-0.500)/0.250)) + 0.65*(1-(0.5+0.5*Foam::tanh((x-0.500)/0.250)));
	//asCells[celli] 	= asCells[celli]*(0.5+0.5*Foam::tanh((x-0.500)/0.250)) + 1*(1-(0.5+0.5*Foam::tanh((x-0.500)/0.250)));
	//muCells[celli] 	= muCells[celli]*(0.5+0.5*Foam::tanh((x-0.500)/0.250)) + 1.131461e-05*(1-(0.5+0.5*Foam::tanh((x-0.500)/0.250)));
	//alphaCells[celli] 	= alphaCells[celli]*(0.5+0.5*Foam::tanh((x-0.500)/0.250)) + 1.9e-6*(1-(0.5+0.5*Foam::tanh((x-0.500)/0.250)));

    }

    //- Boundary conditions
    if (adiabaticMode == true)
    {
        forAll(Z_.boundaryField(), patchi)
        {
            const fvPatchScalarField& pcsi =
                this->Z_.boundaryField()[patchi];
            
            const fvPatchScalarField& pcsiv2 =
                this->Zvar_.boundaryField()[patchi];

            const fvPatchScalarField& pchi_st =
                this->chi_st_.boundaryField()[patchi];

            fvPatchScalarField& pmu =
                this->mu_favre_.boundaryFieldRef()[patchi];

            fvPatchScalarField& palpha =
                this->alpha_favre_.boundaryFieldRef()[patchi];

            fvPatchScalarField& prho =
                this->density_reynolds_.boundaryFieldRef()[patchi];

            fvPatchScalarField& pt = this->T_.boundaryFieldRef()[patchi];

            fvPatchScalarField& pdefectSt = this->defectSt_.boundaryFieldRef()[patchi];             //Minh

            fvPatchScalarField& pas = this->as_.boundaryFieldRef()[patchi];

            fvPatchScalarField& pMW =
                this->MW_.boundaryFieldRef()[patchi];

            fvPatchScalarField& pCp =
                this->Cp_.boundaryFieldRef()[patchi];

            fvPatchScalarField& plambda =
                this->lambda_.boundaryFieldRef()[patchi];


            fvPatchScalarField& pkappa1 =
                this->kappa1_reynolds_.boundaryFieldRef()[patchi];               

            fvPatchScalarField& pkappa2 =
                this->kappa2_reynolds_.boundaryFieldRef()[patchi];              // Rui to have 5 bands WSGG

            fvPatchScalarField& pkappa3 =
                this->kappa3_reynolds_.boundaryFieldRef()[patchi];

            fvPatchScalarField& pkappa4 =
                this->kappa4_reynolds_.boundaryFieldRef()[patchi];


            fvPatchScalarField& pEm1 =
                this->Em1_reynolds_.boundaryFieldRef()[patchi];               //Minh

            fvPatchScalarField& pEm2 =
                this->Em2_reynolds_.boundaryFieldRef()[patchi];

            fvPatchScalarField& pEm3 =
                this->Em3_reynolds_.boundaryFieldRef()[patchi];

            fvPatchScalarField& pEm4 =
                this->Em4_reynolds_.boundaryFieldRef()[patchi];


            fvPatchScalarField& pSumMH =
                this->SumMH_reynolds_.boundaryFieldRef()[patchi];               //Minh

            // fvPatchScalarField& pProdProgVar =
            //     this->ProdProgVar_reynolds_.boundaryFieldRef()[patchi];               //Minh

            fvPatchScalarField& pMolarFracN2 =
                this->MolarFracN2_favre_.boundaryFieldRef()[patchi];               //Rui

            fvPatchScalarField& pMolarFracO2 =
                this->MolarFracO2_favre_.boundaryFieldRef()[patchi];               //Rui

            forAll(pcsi, facei)
            {

                double max_chi = max(small_chi_st, pchi_st[facei]);

                //- Pure oxidizer
                if (pcsi[facei]<=0.)
                {
                    deltaHStExtractionObject.GetMeanValues(0.0, 0.0,max_chi,defect, defectSt); //Minh
                    flamelets_library.GetMeanValues
                    (
                        0.,
                        0.,
                        max_chi,
                        defectSt,                   // Old: defect, Minh: defectSt
                        extracted
                    );
                }

                //- Pure fuel
                else if (pcsi[facei]>=1.)
                {
                    deltaHStExtractionObject.GetMeanValues(1.0, 0.0,max_chi,defect, defectSt); //Minh
                    flamelets_library.GetMeanValues
                    (
                        1.,
                        0.,
                        max_chi,
                        defectSt,                   // Old: defect, Minh: defectSt
                        extracted
                    );
                }

                //- Mixture
                else
                {
                    Zvar_normalized =
                        pcsiv2[facei] / (pcsi[facei]*(1.-pcsi[facei]));

                    if (Zvar_normalized >= 0.98)
                    {
                        deltaHStExtractionObject.GetMeanValues(pcsi[facei], 0.98,max_chi,defect, defectSt); //Minh
                        flamelets_library.GetMeanValues
                        (
                            pcsi[facei],
                            0.98,
                            max_chi,
                            defectSt,               // Old: defect, Minh: defectSt
                            extracted
                        );
                    }
                    else if (Zvar_normalized < 0.)
                    {
                        deltaHStExtractionObject.GetMeanValues(pcsi[facei], 0.0, max_chi,defect, defectSt); //Minh
                        flamelets_library.GetMeanValues
                        (
                            pcsi[facei],
                            0.00,
                            max_chi,
                            defectSt,               // Old: defect, Minh: defectSt
                            extracted
                        );
                    }
                    else
                    {
                        deltaHStExtractionObject.GetMeanValues(pcsi[facei], Zvar_normalized, max_chi,defect, defectSt); //Minh
                        flamelets_library.GetMeanValues
                        (
                            pcsi[facei],
                            Zvar_normalized,
                            max_chi,
                            defectSt,                  // Old: defect, Minh: defectSt
                            extracted
                        );
                    }
                }

                pdefectSt[facei] = defectSt;                // Minh            

                pt[facei] = extracted[1];

                prho[facei] = extracted[2];

                pas[facei] = extracted[3];

		        MWCells[facei] = extracted[4] / 1000.0;  // g/mole to kg/mole

		        CpCells[facei] = extracted[5];

		        lambdaCells[facei] = extracted[6];

                pmu[facei] = extracted[7];

                palpha[facei] = extracted[8];
    
                pkappa1[facei] = extracted[9];          // Minh
                pkappa2[facei] = extracted[10];            // Rui to have 5 bands WSGG
                pkappa3[facei] = extracted[11];
                pkappa4[facei] = extracted[12];           
        
                pEm1[facei] = extracted[13];          //Minh
                pEm2[facei] = extracted[14];
                pEm3[facei] = extracted[15];
                pEm4[facei] = extracted[16];           

                pSumMH[facei] = extracted[17];           // Minh

                // pProdProgVar[facei] = extracted[16];              //Minh

                pMolarFracN2[facei] = extracted[18];              //Rui
                pMolarFracO2[facei] = extracted[19];              //Rui



		// Relaxation until 5 second of properties T, rho, as, mu, alpha
        	//scalar x = this->time().value();
        	//pt[facei]    = pt[facei]*(0.5+0.5*Foam::tanh((x-0.500)/0.250)) + 294.85*(1-(0.5+0.5*Foam::tanh((x-0.500)/0.250)));
       		//prho[facei]  = prho[facei]*(0.5+0.5*Foam::tanh((x-0.500)/0.250)) + 0.65*(1-(0.5+0.5*Foam::tanh((x-0.500)/0.250)));
      		//pas[facei]  = pas[facei]*(0.5+0.5*Foam::tanh((x-0.500)/0.250)) + 1*(1-(0.5+0.5*Foam::tanh((x-0.500)/0.250)));
      		//pmu[facei]  = pmu[facei]*(0.5+0.5*Foam::tanh((x-0.500)/0.250)) + 1.131461e-05*(1-(0.5+0.5*Foam::tanh((x-0.500)/0.250)));
       		//palpha[facei]       = palpha[facei]*(0.5+0.5*Foam::tanh((x-0.500)/0.250)) + 1.9e-6*(1-(0.5+0.5*Foam::tanh((x-0.500)/0.250)));

		

            }
        }
    }
    else
    {
        forAll(Z_.boundaryField(), patchi)
        {

            if (patch_type_T[patchi] == 0)
            {
                const fvPatchScalarField& pcsi =
                    this->Z_.boundaryField()[patchi];

                const fvPatchScalarField& pcsiv2 =
                    this->Zvar_.boundaryField()[patchi];
                
                const fvPatchScalarField& pchi_st =
                    this->chi_st_.boundaryField()[patchi];

                const fvPatchScalarField& ph =
                    this->H_.boundaryField()[patchi];

                fvPatchScalarField& prho =
                    this->density_reynolds_.boundaryFieldRef()[patchi];

                fvPatchScalarField& pmu =
                    this->mu_favre_.boundaryFieldRef()[patchi];

                fvPatchScalarField& palpha =
                    this->alpha_favre_.boundaryFieldRef()[patchi];

                fvPatchScalarField& pas =
                    this->as_.boundaryFieldRef()[patchi];

                fvPatchScalarField& pt = this->T_.boundaryFieldRef()[patchi];

                fvPatchScalarField& pdefectSt = this->defectSt_.boundaryFieldRef()[patchi];     //Minh 


                fvPatchScalarField& pMW =
                	this->MW_.boundaryFieldRef()[patchi];

            	fvPatchScalarField& pCp =
                	this->Cp_.boundaryFieldRef()[patchi];

            	fvPatchScalarField& plambda =
                	this->lambda_.boundaryFieldRef()[patchi];


                fvPatchScalarField& pkappa1 =
                    this->kappa1_reynolds_.boundaryFieldRef()[patchi];                     

                fvPatchScalarField& pkappa2 =
                    this->kappa2_reynolds_.boundaryFieldRef()[patchi];

                fvPatchScalarField& pkappa3 =
                    this->kappa3_reynolds_.boundaryFieldRef()[patchi];

                fvPatchScalarField& pkappa4 =
                    this->kappa4_reynolds_.boundaryFieldRef()[patchi];


                fvPatchScalarField& pEm1 =
                    this->Em1_reynolds_.boundaryFieldRef()[patchi];                  //Minh

                fvPatchScalarField& pEm2 =
                    this->Em2_reynolds_.boundaryFieldRef()[patchi]; 

                fvPatchScalarField& pEm3 =
                    this->Em3_reynolds_.boundaryFieldRef()[patchi]; 

                fvPatchScalarField& pEm4 =
                    this->Em4_reynolds_.boundaryFieldRef()[patchi]; 



                fvPatchScalarField& pSumMH =
                    this->SumMH_reynolds_.boundaryFieldRef()[patchi];                     //Minh

                // fvPatchScalarField& pProdProgVar =
                //     this->ProdProgVar_reynolds_.boundaryFieldRef()[patchi];                  //Minh

                fvPatchScalarField& pMolarFracN2 =
                    this->MolarFracN2_favre_.boundaryFieldRef()[patchi];                  //Rui

                fvPatchScalarField& pMolarFracO2 =
                    this->MolarFracO2_favre_.boundaryFieldRef()[patchi];                  //Rui



                    

                forAll(pcsi, facei)
                {

                    double max_chi = max(small_chi_st, pchi_st[facei]);

                    defect =
                        ph[facei] - (HOxidizer+pcsi[facei]*(HFuel-HOxidizer));

                    //- Pure oxidizer
                    if (pcsi[facei]<=0.)
                    {
                        deltaHStExtractionObject.GetMeanValues(0., 0., max_chi,defect, defectSt); //Minh
                        flamelets_library.GetMeanValues
                        (
                            0.,
                            0.,
                            max_chi,
                            defectSt,               // Old: defect, Minh: defectSt
                            extracted
                        );
                    }

                    //- Pure fuel
                    else if (pcsi[facei]>=1.)
                    {
                        deltaHStExtractionObject.GetMeanValues(1.0, 0., max_chi,defect, defectSt); //Minh
                        flamelets_library.GetMeanValues
                        (
                            1.,
                            0.,
                            max_chi,
                            defectSt,               // Old: defect, Minh: defectSt
                            extracted
                        );
                    }

                    //- Mixture
                    else
                    {
                        Zvar_normalized
                            = pcsiv2[facei] / (pcsi[facei]*(1.-pcsi[facei]));

                        if (Zvar_normalized >= 0.98)
                        {
                            deltaHStExtractionObject.GetMeanValues(pcsi[facei], 0.98, max_chi,defect, defectSt); //Minh
                            flamelets_library.GetMeanValues
                            (
                                pcsi[facei],
                                0.98,
                                max_chi,
                                defectSt,                       // Old: defect, Minh: defectSt
                                extracted
                            );
                        }
                        else if (Zvar_normalized < 0.)
                        {
                            deltaHStExtractionObject.GetMeanValues(pcsi[facei], 0.0, max_chi,defect, defectSt); //Minh
                            flamelets_library.GetMeanValues
                            (
                                pcsi[facei],
                                0.00,
                                max_chi,
                                defectSt,                       // Old: defect, Minh: defectSt
                                extracted
                            );
                        }
                        else
                        {
                            deltaHStExtractionObject.GetMeanValues(pcsi[facei], Zvar_normalized, max_chi,defect, defectSt); //Minh
                            flamelets_library.GetMeanValues
                            (
                                pcsi[facei],
                                Zvar_normalized,
                                max_chi,
                                defectSt,                   // Old: defect, Minh: defectSt
                                extracted
                            );
                        }
                    }

                    pdefectSt [facei] = defectSt;               //Minh

                    pt[facei] = extracted[1];
                    
                    prho[facei] = extracted[2];
                    
                    pas[facei] = extracted[3];

			        MWCells[facei] = extracted[4] / 1000.0;  // g/mole to kg/mole

			        CpCells[facei] = extracted[5];

			        lambdaCells[facei] = extracted[6];

                    pmu[facei] = extracted[7];
                    
                    palpha[facei] = extracted[8];

                    pkappa1[facei] = extracted[9];               //Rui
                    pkappa2[facei] = extracted[10];
                    pkappa3[facei] = extracted[11];
                    pkappa4[facei] = extracted[12];

                    pEm1[facei] = extracted[13];                  //Rui
                    pEm2[facei] = extracted[14];
                    pEm3[facei] = extracted[15];
                    pEm4[facei] = extracted[16];

                    pSumMH[facei] = extracted[17];               //Minh

                    // pProdProgVar[facei] = extracted[16];                  //Minh
                    pMolarFracN2[facei] = extracted[18];                  //Rui
                    pMolarFracO2[facei] = extracted[19];                  //Rui


		    // Relaxation until 5 second of properties T, rho, as, mu, alpha
                //scalar x = this->time().value();
		//pt[facei]    = pt[facei]*(0.5+0.5*Foam::tanh((x-0.500)/0.250)) + 294.75*(1-(0.5+0.5*Foam::tanh((x-0.500)/0.250)));
                //prho[facei]  = prho[facei]*(0.5+0.5*Foam::tanh((x-0.500)/0.250)) + 0.65*(1-(0.5+0.5*Foam::tanh((x-0.500)/0.250)));
                //pas[facei]  = pas[facei]*(0.5+0.5*Foam::tanh((x-0.500)/0.250)) + 1*(1-(0.5+0.5*Foam::tanh((x-0.500)/0.250)));
                //pmu[facei]  = pmu[facei]*(0.5+0.5*Foam::tanh((x-0.500)/0.250)) + 1.131461e-05*(1-(0.5+0.5*Foam::tanh((x-0.500)/0.250)));
                //palpha[facei]       = palpha[facei]*(0.5+0.5*Foam::tanh((x-0.500)/0.250)) + 1.9e-6*(1-(0.5+0.5*Foam::tanh((x-0.500)/0.250)));


                }
            }

            //- Added for friendly handling
            //- fixed temperature and fixed enhalpy :: INLETS
            else if 
            (
                patch_type_T[patchi] == 1
             && patch_type_H[patchi] == 1
             && patch_type_Z[patchi] == 1
            )
            {
                const fvPatchScalarField& pcsi
                    = this->Z_.boundaryField()[patchi];

                const fvPatchScalarField& pcsiv2
                    = this->Zvar_.boundaryField()[patchi];

                const fvPatchScalarField& pchi_st
                    = this->chi_st_.boundaryField()[patchi];

                fvPatchScalarField& prho
                    = this->density_reynolds_.boundaryFieldRef()[patchi];

                fvPatchScalarField& pas
                    = this->as_.boundaryFieldRef()[patchi];

                fvPatchScalarField& pmu
                    = this->mu_favre_.boundaryFieldRef()[patchi];

                fvPatchScalarField& palpha
                    = this->alpha_favre_.boundaryFieldRef()[patchi];

                fvPatchScalarField& pdefectSt = this->defectSt_.boundaryFieldRef()[patchi];     //Minh


                fvPatchScalarField& pMW =
                	this->MW_.boundaryFieldRef()[patchi];

            	fvPatchScalarField& pCp =
                	this->Cp_.boundaryFieldRef()[patchi];

            	fvPatchScalarField& plambda =
                	this->lambda_.boundaryFieldRef()[patchi];


                fvPatchScalarField& pkappa1
                    = this->kappa1_reynolds_.boundaryFieldRef()[patchi];                         

                fvPatchScalarField& pkappa2
                    = this->kappa2_reynolds_.boundaryFieldRef()[patchi];

                fvPatchScalarField& pkappa3
                    = this->kappa3_reynolds_.boundaryFieldRef()[patchi];

                fvPatchScalarField& pkappa4
                    = this->kappa4_reynolds_.boundaryFieldRef()[patchi];



                fvPatchScalarField& pEm1
                    = this->Em1_reynolds_.boundaryFieldRef()[patchi];                            //Minh

                fvPatchScalarField& pEm2
                    = this->Em2_reynolds_.boundaryFieldRef()[patchi];

                fvPatchScalarField& pEm3
                    = this->Em3_reynolds_.boundaryFieldRef()[patchi];

                fvPatchScalarField& pEm4
                    = this->Em4_reynolds_.boundaryFieldRef()[patchi];



                fvPatchScalarField& pSumMH
                    = this->SumMH_reynolds_.boundaryFieldRef()[patchi];                         //Minh

                // fvPatchScalarField& pProdProgVar
                //     = this->ProdProgVar_reynolds_.boundaryFieldRef()[patchi];                            //Minh

                fvPatchScalarField& pMolarFracN2
                    = this->MolarFracN2_favre_.boundaryFieldRef()[patchi];                            //Minh

                fvPatchScalarField& pMolarFracO2
                    = this->MolarFracO2_favre_.boundaryFieldRef()[patchi];                            //Minh


                    

                forAll(pcsi, facei)
                {
                    defect = 0.;

                    double max_chi = max(small_chi_st, pchi_st[facei]);

                    //- Pure oxidizer
                    if (pcsi[facei]<=0.)
                    {
                        deltaHStExtractionObject.GetMeanValues(0.0, 0.0, max_chi,defect, defectSt); //Minh
                        flamelets_library.GetMeanValues
                        (
                            0.,
                            0.,
                            max_chi,
                            defectSt,                       // Old: defect, Minh: defectSt
                            extracted
                        );
                    }

                    //- Pure fuel
                    else if (pcsi[facei]>=1.)
                    {
                        deltaHStExtractionObject.GetMeanValues(1.0, 0.0, max_chi,defect, defectSt); //Minh
                        flamelets_library.GetMeanValues
                        (
                            1.,
                            0.,
                            max_chi,
                            defectSt,                       // Old: defect, Minh: defectSt
                            extracted
                        );
                    }

                    //- Mixture
                    else
                    {
                        Zvar_normalized
                            = pcsiv2[facei] / (pcsi[facei]*(1.-pcsi[facei]));

                        if (Zvar_normalized >= 0.98)
                        {
                            deltaHStExtractionObject.GetMeanValues(pcsi[facei], 0.98, max_chi,defect, defectSt); //Minh
                            flamelets_library.GetMeanValues
                            (
                                pcsi[facei],
                                0.98,
                                max_chi,
                                defectSt,                   // Old: defect, Minh: defectSt
                                extracted
                            );
                        }
                        else if (Zvar_normalized < 0.)
                        {
                            deltaHStExtractionObject.GetMeanValues(pcsi[facei], 0.0, max_chi,defect, defectSt); //Minh
                            flamelets_library.GetMeanValues
                            (
                                pcsi[facei],
                                0.00,
                                max_chi,
                                defectSt,                   // Old: defect, Minh: defectSt
                                extracted
                            );
                        }
                        else
                        {
                            deltaHStExtractionObject.GetMeanValues(pcsi[facei], Zvar_normalized, max_chi,defect, defectSt); //Minh
                            flamelets_library.GetMeanValues
                            (
                                pcsi[facei],
                                Zvar_normalized,
                                max_chi,
                                defectSt,                   // Old: defect, Minh: defectSt
                                extracted
                            );
                        }
                    }

                    pdefectSt[facei] = defectSt;            //Minh

                    prho[facei] = extracted[2];

                    pas[facei] = extracted[3];

			        MWCells[facei] = extracted[4] / 1000.0;  // g/mole to kg/mole

			        CpCells[facei] = extracted[5];

			        lambdaCells[facei] = extracted[6];

                    pmu[facei] = extracted[7];

                    palpha[facei] = extracted[8];

                    pkappa1[facei] = extracted[9];           //Minh
                    pkappa2[facei] = extracted[10];
                    pkappa3[facei] = extracted[11];
                    pkappa4[facei] = extracted[12];

                    pEm1[facei] = extracted[13];              //Minh
                    pEm2[facei] = extracted[14];
                    pEm3[facei] = extracted[15];
                    pEm4[facei] = extracted[16];


                    pSumMH[facei] = extracted[17];               //Minh

                    // pProdProgVar[facei] = extracted[16];                  //Minh

                    pMolarFracN2[facei] = extracted[18];               //Rui
                    pMolarFracO2[facei] = extracted[19];               //Rui


 		    // Relaxation until 5 second of properties T, rho, as, mu, alpha
                //scalar x = this->time().value();
                //prho[facei]  = prho[facei]*(0.5+0.5*Foam::tanh((x-0.500)/0.250)) + 0.65*(1-(0.5+0.5*Foam::tanh((x-0.500)/0.250)));
                //pas[facei]  = pas[facei]*(0.5+0.5*Foam::tanh((x-0.500)/0.250)) + 1*(1-(0.5+0.5*Foam::tanh((x-0.500)/0.250)));
                //pmu[facei]  = pmu[facei]*(0.5+0.5*Foam::tanh((x-0.500)/0.250)) + 1.131461e-05*(1-(0.5+0.5*Foam::tanh((x-0.500)/0.250)));
                //palpha[facei]       = palpha[facei]*(0.5+0.5*Foam::tanh((x-0.500)/0.250)) + 1.9e-6*(1-(0.5+0.5*Foam::tanh((x-0.500)/0.250)));



		
                }
            }

            //- Added for enthalpydefect due to fixedTemperature on walls
            //- fixed temperature and fixedEnthalpie :: no inlet
            else if
            (
                patch_type_T[patchi] == 1
             && patch_type_H[patchi] == 1
             && patch_type_Z[patchi] == 0
            )
            {

                const fvPatchScalarField& pcsi
                    = this->Z_.boundaryField()[patchi];

                const fvPatchScalarField& pcsiv2
                    = this->Zvar_.boundaryField()[patchi];

                const fvPatchScalarField& pchi_st
                    = this->chi_st_.boundaryField()[patchi];

                fvPatchScalarField& ph
                    = this->H_.boundaryFieldRef()[patchi];

                fvPatchScalarField& pt
                    = this->T_.boundaryFieldRef()[patchi];

                fvPatchScalarField& prho
                    = this->density_reynolds_.boundaryFieldRef()[patchi];

                fvPatchScalarField& pas
                    = this->as_.boundaryFieldRef()[patchi];

                fvPatchScalarField& pmu
                    = this->mu_favre_.boundaryFieldRef()[patchi];

                fvPatchScalarField& palpha
                    = this->alpha_favre_.boundaryFieldRef()[patchi];

                fvPatchScalarField& pdefectSt = this->defectSt_.boundaryFieldRef()[patchi];     //Minh


                fvPatchScalarField& pMW =
                	this->MW_.boundaryFieldRef()[patchi];

            	fvPatchScalarField& pCp =
                	this->Cp_.boundaryFieldRef()[patchi];

            	fvPatchScalarField& plambda =
                	this->lambda_.boundaryFieldRef()[patchi];



                fvPatchScalarField& pkappa1
                    = this->kappa1_reynolds_.boundaryFieldRef()[patchi];                     

                fvPatchScalarField& pkappa2
                    = this->kappa2_reynolds_.boundaryFieldRef()[patchi];

                fvPatchScalarField& pkappa3
                    = this->kappa3_reynolds_.boundaryFieldRef()[patchi];

                fvPatchScalarField& pkappa4
                    = this->kappa4_reynolds_.boundaryFieldRef()[patchi];


                fvPatchScalarField& pEm1
                    = this->Em1_reynolds_.boundaryFieldRef()[patchi];                    

                fvPatchScalarField& pEm2
                    = this->Em2_reynolds_.boundaryFieldRef()[patchi]; 

                fvPatchScalarField& pEm3
                    = this->Em3_reynolds_.boundaryFieldRef()[patchi]; 

                fvPatchScalarField& pEm4
                    = this->Em4_reynolds_.boundaryFieldRef()[patchi]; 



                fvPatchScalarField& pSumMH
                    = this->SumMH_reynolds_.boundaryFieldRef()[patchi];                     //Minh


                // fvPatchScalarField& pProdProgVar
                //     = this->ProdProgVar_reynolds_.boundaryFieldRef()[patchi];                    //Minh

                fvPatchScalarField& pMolarFracN2
                    = this->MolarFracN2_favre_.boundaryFieldRef()[patchi];                    //Rui

                fvPatchScalarField& pMolarFracO2
                    = this->MolarFracO2_favre_.boundaryFieldRef()[patchi];                    //Rui




                    

                forAll(pcsi, facei)
                {

                    scalar max_chi = max(small_chi_st, pchi_st[facei]);

                    //- Pure oxidizer
                    if (pcsi[facei]<=0.)
                    {
                        defect =
                            flamelets_library.GetEnthalpyDefectFromTemperature
                                (
                                    0.,
                                    0.,
                                    max_chi,
                                    pt[facei]
                                );

                        ph[facei] = defect + HOxidizer;
                        deltaHStExtractionObject.GetMeanValues(0., 0., max_chi,defect, defectSt); //Minh
                        flamelets_library.GetMeanValues
                        (
                            0.,
                            0.,
                            max_chi,
                            defectSt,                       // Old: defect, Minh: defectSt
                            extracted
                        );
                    }

                    //- Pure fuel
                    else if (pcsi[facei]>=1.)
                    {
                        defect = 
                            flamelets_library.GetEnthalpyDefectFromTemperature
                            (
                                1.,
                                0.,
                                max_chi,
                                pt[facei]
                            );

                        ph[facei] = defect + HFuel;
                            
                        deltaHStExtractionObject.GetMeanValues(1., 0., max_chi,defect, defectSt); //Minh
                        flamelets_library.GetMeanValues
                        (
                            1.,
                            0.,
                            max_chi,
                            defectSt,                        // Old: defect, Minh: defectSt
                            extracted
                        );
                    }

                    //- Mixture
                    else
                    {
                        Zvar_normalized =
                            pcsiv2[facei] / (pcsi[facei]*(1.-pcsi[facei]));

                        if (Zvar_normalized >= 0.98)
                        {
                            defect =
                                flamelets_library.GetEnthalpyDefectFromTemperature
                                (
                                    pcsi[facei],
                                    0.98,
                                    max_chi,
                                    pt[facei]
                                );

                            ph[facei] =
                                defect
                              + (HOxidizer+pcsi[facei]*(HFuel-HOxidizer));
                    
                            deltaHStExtractionObject.GetMeanValues(pcsi[facei], 0.98, max_chi,defect, defectSt); //Minh
                            flamelets_library.GetMeanValues
                            (
                                pcsi[facei],
                                0.98,
                                max_chi,                  
                                defectSt,                   // Old: defect, Minh: defectSt
                                extracted
                            );
                        }
                        else if (Zvar_normalized < 0.)
                        {
                            defect =
                                flamelets_library.GetEnthalpyDefectFromTemperature
                                (
                                    pcsi[facei],
                                    0.00,
                                    max_chi,
                                    pt[facei]
                                );

                            ph[facei] = 
                                defect
                              + (HOxidizer+pcsi[facei]*(HFuel-HOxidizer));
        
                            deltaHStExtractionObject.GetMeanValues(pcsi[facei], 0.0, max_chi,defect, defectSt); //Minh
                            
                            flamelets_library.GetMeanValues
                            (
                                pcsi[facei],
                                0.00,
                                max_chi,
                                defectSt,                   // Old: defect, Minh: defectSt
                                extracted
                            );
                        }
                        else
                        {
                            defect =
                                flamelets_library.GetEnthalpyDefectFromTemperature
                                (
                                    pcsi[facei],
                                    Zvar_normalized,
                                    max_chi,
                                    pt[facei]
                                );

                            ph[facei] =
                                defect
                              + (HOxidizer+pcsi[facei]*(HFuel-HOxidizer));
            
                            deltaHStExtractionObject.GetMeanValues(pcsi[facei], Zvar_normalized, max_chi,defect, defectSt); //Minh
                            flamelets_library.GetMeanValues
                            (
                                pcsi[facei],
                                Zvar_normalized,
                                max_chi,
                                defectSt,                   // Old: defect, Minh: defectSt
                                extracted
                            );
                        }
                    }

                    pdefectSt[facei] = defectSt;              //Minh

                    prho[facei] = extracted[2];

                    pas[facei] = extracted[3];

			        MWCells[facei] = extracted[4] / 1000.0;  // g/mole to kg/mole

			        CpCells[facei] = extracted[5];

			        lambdaCells[facei] = extracted[6];

                    pmu[facei] = extracted[7];

                    palpha[facei] = extracted[8];

                    pkappa1[facei] = extracted[9];               //Rui
                    pkappa2[facei] = extracted[10];
                    pkappa3[facei] = extracted[11];
                    pkappa4[facei] = extracted[12];

                    pEm1[facei] = extracted[13];                  //Rui
                    pEm2[facei] = extracted[14];
                    pEm3[facei] = extracted[15];
                    pEm4[facei] = extracted[16];
 

                    pSumMH[facei] = extracted[17];               //Minh

                    // pProdProgVar[facei] = extracted[16];                  //Minh

                    pMolarFracN2[facei] = extracted[18];               //Rui
                    pMolarFracO2[facei] = extracted[19];               //Rui




 		    // Relaxation until 5 second of properties T, rho, as, mu, alpha
                //scalar x = this->time().value();
                //prho[facei]  = prho[facei]*(0.5+0.5*Foam::tanh((x-0.500)/0.250)) + 0.65*(1-(0.5+0.5*Foam::tanh((x-0.500)/0.250)));
                //pas[facei]  = pas[facei]*(0.5+0.5*Foam::tanh((x-0.500)/0.250)) + 1*(1-(0.5+0.5*Foam::tanh((x-0.500)/0.250)));
                //pmu[facei]  = pmu[facei]*(0.5+0.5*Foam::tanh((x-0.500)/0.250)) + 1.131461e-05*(1-(0.5+0.5*Foam::tanh((x-0.500)/0.250)));
                //palpha[facei]       = palpha[facei]*(0.5+0.5*Foam::tanh((x-0.500)/0.250)) + 1.9e-6*(1-(0.5+0.5*Foam::tanh((x-0.500)/0.250)));


                }
            }

            //- Not implemented boundary condition :: give an error
            else
            {
                FatalErrorIn
                (
                    "pdfFlameletPsiReactionThermo<BasicFlameletThermo,"
                    "MixtureType>::update()"
                )
                << "Boundary conditions are wrong: "
                << "fixed temperature BC must be fixed "
                << "enthaplie BC with value 0;"
                << abort(FatalError);
            }
        }
    }
}

template<class BasicFlameletThermo, class MixtureType>
void Foam::pdfFlameletPsiReactionThermo<BasicFlameletThermo, MixtureType>::
updateMassFractions()
{
    std::vector<scalar> extracted(flamelets_library.number_of_species()+1);

    double Zvar_normalized = 0.;

    double defect = 0.;
      
    double defectSt = 0.0;

    const scalarField& Z = this->Z_.internalField();

    const scalarField& Zvar = this->Zvar_.internalField();

    const scalarField& chi_st = this->chi_st_.internalField();

    const scalarField& HCells = this->H_.internalField();

    double small_eps = 1.e-6;

    double small_chi_st = 1.e-8;

    Foam::PtrList<Foam::volScalarField>& Y = this->Y_; // Minh

    volScalarField Yt(0.0*Y[0]);	

    forAll(Z, celli)
    {
        double max_chi = max(small_chi_st,chi_st[celli]);

        if (adiabaticMode == false)
        {
            defect = HCells[celli] - (HOxidizer+Z[celli]*(HFuel-HOxidizer));
        }

        //- Pure oxidizer
        if (Z[celli]<=small_eps)
        {
            deltaHStExtractionObject.GetMeanValues(0., 0., max_chi,defect, defectSt); //Minh
            flamelets_library.ExtractMeanValues
            (
                0.,
                0.,
                max_chi,
                defectSt,                       // Old: defect, Minh: defectSt
                extracted
            );
        }

        //- Pure fuel
        else if (Z[celli]>=(1.-small_eps))
        {
            deltaHStExtractionObject.GetMeanValues(1., 0., max_chi,defect, defectSt); //Minh
            flamelets_library.ExtractMeanValues
            (
                1.,
                0.,
                max_chi,
                defectSt,                       // Old: defect, Minh: defectSt
                extracted
            );
        }

        //- Mixture
        else
        {
            Zvar_normalized = Zvar[celli] / (Z[celli]*(1.-Z[celli]));

            if (Zvar_normalized >= 0.98)
            {
                deltaHStExtractionObject.GetMeanValues(Z[celli], 0.98, max_chi,defect, defectSt); //Minh
                flamelets_library.ExtractMeanValues
                (
                    Z[celli],
                    0.98,
                    max_chi,
                    defectSt,               // Old: defect, Minh: defectSt
                    extracted
                );
            }
            else if (Zvar_normalized < 0.)
            {
                deltaHStExtractionObject.GetMeanValues(Z[celli], 0.00, max_chi,defect, defectSt); //Minh
                flamelets_library.ExtractMeanValues
                (
                    Z[celli],
                    0.00,
                    max_chi,
                    defectSt,                   // Old: defect, Minh: defectSt
                    extracted
                );
            }
            else
            {
                deltaHStExtractionObject.GetMeanValues(Z[celli], Zvar_normalized, max_chi,defect, defectSt); //Minh
                flamelets_library.ExtractMeanValues
                (
                    Z[celli],
                    Zvar_normalized,
                    max_chi,
                    defectSt,                   // Old: defect, Minh: defectSt
                    extracted
                );
            }
        }

        for(int j=0;j<flamelets_library.number_of_species()+1;j++)
        {
            if(j<(flamelets_library.number_of_species() -1))
            {
                omega_[j].ref()[celli] = extracted[j+1];
		Y[j].ref()[celli] = extracted[j+1];		//Minh
		Yt.ref()[celli] += Y[j].ref()[celli];
            }
	
	    if(j == (flamelets_library.number_of_species()-1))
	    {
		Y[flamelets_library.number_of_species()-1].ref()[celli] = scalar(1.0)- Yt.ref()[celli];
            }	 	
        }
    }

    forAll(Z_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pZ = this->Z_.boundaryField()[patchi];
        
        const fvPatchScalarField& pZvar = this->Zvar_.boundaryField()[patchi];

        const fvPatchScalarField& ph = this->H_.boundaryField()[patchi];

        const fvPatchScalarField& pchi_st =
            this->chi_st_.boundaryField()[patchi];

        forAll(pZ, facei)
        {

            double max_chi = max(small_chi_st, pchi_st[facei]);

            if (adiabaticMode == false)
                defect = ph[facei] - (HOxidizer+pZ[facei]*(HFuel-HOxidizer));

            //- Pure oxidizer
            if (pZ[facei]<=small_eps)
            {
                deltaHStExtractionObject.GetMeanValues(0.0, 0.0, max_chi,defect, defectSt); //Minh
                flamelets_library.ExtractMeanValues
                (
                    0.,
                    0., 
                    max_chi, 
                    defectSt,                   // Old: defect, Minh: defectSt 
                    extracted
                );
            }

            //- Pure fuel
            else if (pZ[facei]>=(1.-small_eps))
            {
                deltaHStExtractionObject.GetMeanValues(1.0, 0.0, max_chi,defect, defectSt); //Minh
                flamelets_library.ExtractMeanValues
                (
                    1., 
                    0., 
                    max_chi, 
                    defectSt,                           // Old: defect, Minh: defectSt
                    extracted
                );
            }

            //- Mixture
            else
            {
                Zvar_normalized = pZvar[facei] / (pZ[facei]*(1.-pZ[facei]));

                if (Zvar_normalized >= 0.98)
                {
                    deltaHStExtractionObject.GetMeanValues(pZ[facei], 0.98, max_chi,defect, defectSt); //Minh
                    flamelets_library.ExtractMeanValues
                    (
                        pZ[facei], 
                        0.98, 
                        max_chi, 
                        defectSt,                       // Old: defect, Minh: defectSt 
                        extracted
                    );
                }
                else if (Zvar_normalized < 0.)
                {
                    deltaHStExtractionObject.GetMeanValues(pZ[facei], 0.0, max_chi,defect, defectSt); //Min
                    flamelets_library.ExtractMeanValues
                    (
                        pZ[facei], 
                        0.00, 
                        max_chi, 
                        defectSt,                           // Old: defect, Minh: defectSt 
                        extracted
                   );
                }
                else
                {
                    deltaHStExtractionObject.GetMeanValues(pZ[facei], Zvar_normalized, max_chi,defect, defectSt); //Minh
                    flamelets_library.ExtractMeanValues
                    (
                        pZ[facei], 
                        Zvar_normalized, 
                        max_chi,
                        defectSt,                               // Old: defect, Minh: defectSt 
                        extracted
                    );
                }
            }

            for(int j=0;j<flamelets_library.number_of_species()+1;j++)
            {
                if(j<(flamelets_library.number_of_species()-1))
                {
                    omega_[j].boundaryFieldRef()[patchi][facei] =
                        extracted[j+1];
		    Y[j].boundaryFieldRef()[patchi][facei] =
                        extracted[j+1];		//Minh
		    Yt.boundaryFieldRef()[patchi][facei] += Y[j].boundaryFieldRef()[patchi][facei];
                }

 		if(j==(flamelets_library.number_of_species()-1))
		{
		    Y[flamelets_library.number_of_species()-1].boundaryFieldRef()[patchi][facei] = scalar(1.00) - Yt.boundaryFieldRef()[patchi][facei];

		}

            }
        }
    }
}



template<class BasicFlameletThermo, class MixtureType>
void Foam::pdfFlameletPsiReactionThermo<BasicFlameletThermo, MixtureType>::
errorMessage(const string message)
{
    Info << "Class: pdfFlameletPsiReactionThermo" << endl;
    Info << "Error: " << message << endl;
    getchar();
}

template<class BasicFlameletThermo, class MixtureType>
void Foam::pdfFlameletPsiReactionThermo<BasicFlameletThermo, MixtureType>::
infoMessage() const
{
    Info << "/*-------------------------------------*\\\n"
         << "|    Flamelet Thermo initialization     |\n"
         << "|                                       |\n"
         << "|   Rebuild by Tobias Holzmann M.Eng.   |\n"
         << "|          www.Holzmann-cfd.de          |\n"
         << "\\*-------------------------------------*/\n"
         << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::pdfFlameletPsiReactionThermo<BasicPsiThermo, MixtureType>::pdfFlameletPsiReactionThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    heThermo<BasicPsiThermo, MixtureType>(mesh, phaseName),

    Z_
    (
        IOobject
        (
            "Z",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    Zvar_
    (
        IOobject
        (
            "Zvar",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
        // dimensionedScalar("Zvar", dimless,0.0)
    ),

    chi_st_
    (
        IOobject
        (
            "chi_st",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
        // dimensionedScalar("chi_st",dimensionSet(0,0,-1,0,0,0,0) , 0.0)
    ),

    H_
    (
        IOobject
        (
            "H",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    defect_
    (
        IOobject
        (
            "defect",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("defect",dimensionSet(0,2,-2,0,0,0,0) , 0.0)
    ),

    defectSt_
    (
        IOobject
        (
            "defectSt",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("defectSt",dimensionSet(0,2,-2,0,0,0,0) , 0.0)
    ),

    as_
    (
        IOobject
        (
            "as",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("as",dimensionSet(0,-1,0,0,0,0,0) , 0.0)
    ),

    MW_     
    (
        IOobject
        (
            "MW",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("MW",dimensionSet(1,0,0,0,-1,0,0) , 0.0)     // kg/mol, original is g/mol
    ),

    Cp_
    (
        IOobject
        (
            "Cp",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Cp",dimensionSet(0,2,-2,-1,0,0,0) , 0.0)     // m2/s2/K, original is J/kg/K
    ),

    lambda_
    (
        IOobject
        (
            "lambda",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("lambda",dimensionSet(1,1,-3,-1,0,0,0) , 0.0)     // kg*m/s3/K, original is J/s/m/K
    ),

    kappa1_reynolds_
    (
        IOobject
        (
            "kappa1TRI",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("kappa1TRI",dimensionSet(0,-1,0,0,0,0,0) , 0.0)
    ),

    kappa2_reynolds_
    (
        IOobject
        (
            "kappa2TRI",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("kappa2TRI",dimensionSet(0,-1,0,0,0,0,0) , 0.0)
    ),

    kappa3_reynolds_
    (
        IOobject
        (
            "kappa3TRI",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("kappa3TRI",dimensionSet(0,-1,0,0,0,0,0) , 0.0)
    ),

    kappa4_reynolds_
    (
        IOobject
        (
            "kappa4TRI",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("kappa4TRI",dimensionSet(0,-1,0,0,0,0,0) , 0.0)
    ),


    Em1_reynolds_
    (
        IOobject
        (
            "Em1TRI",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("Em1TRI",dimensionSet(1,-1,-3,0,0,0,0) , 0.0)
    ),

    Em2_reynolds_
    (
        IOobject
        (
            "Em2TRI",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("Em2TRI",dimensionSet(1,-1,-3,0,0,0,0) , 0.0)
    ),

    Em3_reynolds_
    (
        IOobject
        (
            "Em3TRI",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("Em3TRI",dimensionSet(1,-1,-3,0,0,0,0) , 0.0)
    ),

    Em4_reynolds_
    (
        IOobject
        (
            "Em4TRI",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("Em4TRI",dimensionSet(1,-1,-3,0,0,0,0) , 0.0)
    ),



    SumMH_reynolds_
    (
        IOobject
        (
            "SumMHFl",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("SumMHFl",dimensionSet(1,-1,-3,0,0,0,0) , 0.0)
    ),

    // ProdProgVar_reynolds_
    // (
    //     IOobject
    //     (
    //         "ProdProgVarFl",
    //         mesh.time().timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh,
    //     dimensionedScalar("ProdProgVarFl",dimensionSet(1,-3,-1,0,0,0,0) , 0.0)
    // ),

    MolarFracN2_favre_
    (
        IOobject
        (
            "MolarFracN2Fl",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("MolarFracN2Fl",dimless, 0.0)
    ),

    MolarFracO2_favre_
    (
        IOobject
        (
            "MolarFracO2Fl",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("MolarFracO2Fl",dimless, 0.0)
    ),

    

    density_reynolds_
    (
        IOobject
        (
            "rho_reynolds",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("density_reynolds",dimensionSet(1,-3,0,0,0,0,0) , 0.0)
    ),

    mu_favre_
    (
        IOobject
        (
            "mu_lam",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("mu_lam",dimensionSet(1,-1,-1,0,0,0,0) , 0.0)
    ),

    alpha_favre_
    (
        IOobject
        (
            "alpha_lam",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("alpha_lam",dimensionSet(0,2,-1,0,0,0,0) , 0.0)
    ),

    adiabaticMode(false),							// initalization false value for the variables
    showFlamelet(false),
    showFlameletLibrary(false)
{
    Info<< "pdfFlameletPsiReactionThermo class object is constructed\n" << endl; 		// Print for checking this constructor is called?
	

    //- Info message
    infoMessage();

    //- Get fixedValue enthalpy boundarys
    Info << "Enthalpy:" << endl;
    forAll(this->H_.boundaryField(), patchi)					// patchi is a bien chay
    {
        if (isA<fixedValueFvPatchScalarField>(this->H_.boundaryField()[patchi]))
        {
            Info<< "     + " << mesh.boundary()[patchi].name()
                << " <fixedValue>" << endl;
            patch_type_H.append(1);
        }
        else
        {
            Info<< "     + " << mesh.boundary()[patchi].name()
                << "" << endl;
            patch_type_H.append(0);
        }
    }

    //- Get fixedValue temperatur boundarys
    Info << endl << "Temperature: " << endl;
    forAll(this->T_.boundaryField(), patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(this->T_.boundaryField()[patchi]))
        {
            Info<< "     + " << mesh.boundary()[patchi].name()
                << " <fixedValue>" << endl;
            patch_type_T.append(1);
        }
        else
        {
            Info << "     + " << mesh.boundary()[patchi].name() << "" << endl;
            patch_type_T.append(0);
        }
    }

    //- Get fuel and oxidizer inlets
    Info << endl << "Mixture fraction:" << endl;
    forAll(this->Z_.boundaryField(), patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(this->Z_.boundaryField()[patchi]))
        {
            Info << "     + " << mesh.boundary()[patchi].name()
                << " <fixedValue>" << endl;
            patch_type_Z.append(1);	
        }
        else
        {
            Info << "     + " << mesh.boundary()[patchi].name() << "" << endl;
            patch_type_Z.append(0);
        }
    }
    Info << endl;

    //- IOFlamelet properties
    IOdictionary flameletsProperties_
    (
        IOobject
        (
            "flameletsProperties",
            Z_.time().constant(),
            Z_.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    Info << "Print flameletsProperties_ object of IOdictionary" << flameletsProperties_ << endl; 	// print flameletsProperties_
    //- Get flamelet properties and set fields and variables
    {
        Switch adiabaticMode_(flameletsProperties_.lookup("adiabaticMode"));
	Info << "Print adiabaticMode_ variable:" << adiabaticMode_ << endl;				// print adiabaticMode_
        Switch showFlamelet_(flameletsProperties_.lookup("showFlamelet"));
        Switch showFlameletLibrary_
        (
            flameletsProperties_.lookup("showFlameletLibrary")
        );
        propertyUpdate =
            readLabel(flameletsProperties_.lookup("propertyUpdate"));

        massFractionsUpdate =
            readLabel(flameletsProperties_.lookup("massFractionsUpdate"));

        adiabaticMode = adiabaticMode_;
        showFlamelet = showFlamelet_;
        showFlameletLibrary = showFlameletLibrary_;
        counter = propertyUpdate;
        counter_mass_fractions = massFractionsUpdate;

        string libraryPath = flameletsProperties_.lookup("libraryPath");
        string chiPDF = flameletsProperties_.lookup("pdf");
        string list_of_species = flameletsProperties_.lookup("species");

        flamelets_library.SetLibraryPath(libraryPath);			//OpenSMOKE_PDF_NonAdiabaticFlamelet_Library flamelets_library;
        flamelets_library.SetSpeciesToExtract(list_of_species);
	

	//- Set the dissipation mode

            //- Set adiabatic mode
            if (adiabaticMode == true)						// Maybe OpenSMOKE_PDF_NonAdiabaticFlamelet_Library is extention OpenSMOKE_PDF_AdiabaticFlamelet_Library class
            {
                Info << "Flamelet thermo mode is <adiabatic>" << endl;
                flamelets_library.SetAdiabaticMode();
            }
            else
            {
                Info << "Flamelet thermo mode is <non-adiabatic>" << endl;
            }

            //- Scalar dissipation rate distribution
            if (chiPDF == "logNormal")
            {
                Info<< "Flamelet thermo use the <log-normal distribution> "
                    << "for the scalar dissipation rate" << endl;

                scalar chi_sigma
                (
                    readScalar(flameletsProperties_.lookup("sigma"))
                );

                label chi_points
                (
                    readScalar(flameletsProperties_.lookup("points"))
                );

                flamelets_library.UnsetExcludeColdFlamelets();
                flamelets_library.SetLogNormalChiDistribution
                (
                    chi_sigma,
                    chi_points
                );
            }
            else if (chiPDF == "dirac")
            {
                Info<< "Flamelet thermo uses the <delta-dirac distribution> "
                    << "for the scalar dissipation rate" << endl;
            }

            //- Set show flamelet mode
            if (showFlamelet == true)
            {
                Info<< "Flamelet thermo shows all single flamelet properties"
                    << endl;

                flamelets_library.SetShowFlamelet();
            }
            else
            {
                Info<< "Flamelet thermo does not show the single "
                    << "flamelet properties" << endl;
            }
            //- Set show flamelet library mode
            if (showFlameletLibrary == true)
            {
                Info<< "Flamelet thermo shows all flamelet library properties"
                    << endl;

                flamelets_library.SetShowFlameletLibrary();
            }
            else
            {
                Info<< "Flamelet thermo does not show the flamelet "
                    << "library properties" << endl;
            }
        //- Initialise the wanted mass fraction fields
        omega_.setSize(flamelets_library.number_of_species()+1);
	Info << "Number of species: " << flamelets_library.number_of_species() << endl; 		// Number of species
        for (int j=0;j<flamelets_library.number_of_species()+1;j++)
        {
            if(j < flamelets_library.number_of_species())
            {
                word name_of_species =
                    "omega_" + flamelets_library.species()[j+1];
		Info << "Name of species" << name_of_species << endl;					// Print name of species
                omega_.set
                (
                    j,
                    new volScalarField
                    (
                        IOobject
                        (
                            name_of_species,
                            mesh.time().timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE
                        ),
                        mesh,
                        dimensionedScalar
                        (
                            "zero",
                            dimensionSet(0,0,0,0,0,0,0),
                            scalar(0.0)
                        )
                    )
                );
	    }
        }
    }

    //- Read all flamelets
    flamelets_library.Read();
    flamelets_library.Summary();

    //- Read the 1st LUT to extract deltaHSt
    deltaHStExtractionObject.Read();

    //- Set the adiabat enthalpy
    Info<< "Thermo flamelet set the adiabat enthalpy for enthalpy "
        << "defect calculation:\n" << endl;
    {
        HOxidizer = flamelets_library.enthalpy_f_oxidizer();
        HFuel = flamelets_library.enthalpy_f_fuel();

        Info<< "     + Adiabat enthalpy fuel: " << HFuel << endl;
        Info<< "     + Adiabat enthalpy oxidizer: " << HOxidizer << "\n"
            << endl;
    }

    //- Extract all variables
    update();

    //- Calculate first time
    calculate();

    // Switch on saving old time
    this->psi_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::pdfFlameletPsiReactionThermo<BasicPsiThermo, MixtureType>::~pdfFlameletPsiReactionThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::pdfFlameletPsiReactionThermo<BasicPsiThermo, MixtureType>::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    //- Force saving of the old-time values
    this->psi_.oldTime();

    if
    (
        counter == propertyUpdate
     || counter_mass_fractions == massFractionsUpdate
    )
    {
        Info<< "Flamelet thermo: \n";
    }

    if (counter == propertyUpdate)
    {
        Info<< "    + Flamelet thermo updates all thermo variables "
            << "using the Look-Up-Table\n";

        update();

        counter = 0;
    }

    if (counter_mass_fractions == massFractionsUpdate)
    {
        Info<< "    + Flamelet thermo updates all mass fractions "
            << " using the Look-Up-Table\n";

        updateMassFractions();

        counter_mass_fractions = 0;
    }
    else if (counter_mass_fractions != massFractionsUpdate && counter == 0)
    {
        Info<< "\n";
    }
        
    calculate();

    ++counter;
    ++counter_mass_fractions;

    if (debug)
    {
        Info<< "    Finished" << endl;
    }
}


// ************************************************************************* //
