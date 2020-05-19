//-  Solving equation of mixture fraction
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

    fvScalarMatrix ZEqn
    (
        fvm::ddt(rho, Z)
      + fvm::div(phi, Z)
      - fvm::laplacian(dEff, Z)
    );
    ZEqn.relax();
    ZEqn.solve(mesh.solver("Z"));
    Z.correctBoundaryConditions();
}

// {
//     //- Solving equation of variance of mixture fraction
//     fvScalarMatrix ZvarEqn
//     (
//         (
//         fvm::ddt(rho, Zvar)
//           + fvm::div(phi, Zvar)
//           - fvm::laplacian(turbulence->mut()/0.5 + thermo.alpha(), Zvar)
//           - 2.0*(turbulence->mut()/0.5)*magSqr(fvc::grad(Z))
//           + 2.0*rho*(turbulence->epsilon()/turbulence->k())*Zvar
//         )
//     );

//     ZvarEqn.relax();

//     ZvarEqn.solve(mesh.solver("Zvar"));
    
//     Zvar.correctBoundaryConditions();
    
//     bound(Zvar, 0.0);
// }

/* if (varianceEquation == false)
{
    //- Solving equation of variance of mixture fraction
    fvScalarMatrix ZvarEqn
    (
        (
	    fvm::ddt(rho, Zvar)
          + fvm::div(phi, Zvar)
          - fvm::laplacian(turbulence->mut()/sigmat, Zvar)
          - Cg*turbulence->mut()*magSqr(fvc::grad(Z))
          + Cd*rho*turbulence->epsilon()/(turbulence->k()+k_small)*Zvar
        )
    );

    ZvarEqn.relax();

    ZvarEqn.solve(mesh.solver("Zvar"));
}
else
{
    Info << "Algebraic equation for Zvar" << endl;

    scalarField& ZvarCells = Zvar.ref();

    ZvarCells =
        Cg*turbulence->mut()*magSqr(fvc::grad(Z))*turbulence->k()
      / (Cd*rho*(turbulence->epsilon()+epsilon_small));
}
 */


//- Updating variance of mixture fraction
scalarField& ZvarCells = Zvar.ref();
// ZvarCells = C_Z*magSqr(sgs.delta())*magSqr(fvc::grad(Z));


ZgradCells = mag(fvc::grad(Z));
deltaCells = mag(sgs.delta());


forAll(ZvarCells, celli)
{   

    // case(1)
    if ( ZgradCells[celli] == 0.0 || Z[celli] == 0.0 || Z[celli] == 1.0 ){
        ZvarCells[celli] = 0.0 ;
    }
    // case(2)
    else if ( ZgradCells[celli] <= 2.0*Z[celli]/deltaCells[celli] && ZgradCells[celli] <= 2.0*(1.0-Z[celli])/deltaCells[celli] )  {
        ZvarCells[celli] = Foam::pow(ZgradCells[celli], 2.0) * Foam::pow(deltaCells[celli], 2.0) / 12.0 ;
    }
    // case(3a)
    else if ( ZgradCells[celli] >= 2.0*Z[celli]/deltaCells[celli] && ZgradCells[celli] <= 1.0/deltaCells[celli]/(2.0*Z[celli]) ) {
        ZvarCells[celli] = 2.0/3.0*Foam::pow(Z[celli], 1.5)*Foam::pow(2.0*ZgradCells[celli]*deltaCells[celli],0.5) - Foam::pow(Z[celli],2.0) ;
    }

    // case (3b)
    else if ( ZgradCells[celli] >= 2.0*(1.0-Z[celli])/deltaCells[celli] && ZgradCells[celli] <= 1.0/(deltaCells[celli]*2.0*(1.0-Z[celli])) ) {
        ZvarCells[celli] = 2.0/3.0*Foam::pow((1.0-Z[celli]),1.5) * Foam::pow(2.0*ZgradCells[celli]*deltaCells[celli],0.5) - Foam::pow((1-Z[celli]),2.0);
    }

    // // case (4)
    // else if ( Z[celli]<=0.5 && ZgradCells[celli] >= 1.0/(deltaCells[celli]*(2.0*Z[celli])) ) {
    //     ZvarCells[celli] = Z[celli]*(1.0-Z[celli])-1.0/(6.0*ZgradCells[celli]*deltaCells[celli]) ;
    // }

    // // case (4)
    // else if (Z[celli] >=0.5 && ZgradCells[celli] >= 1.0/(deltaCells[celli]*2.0*(1.0-Z[celli])) ) {
    //     ZvarCells[celli] = Z[celli]*(1.0-Z[celli])-1.0/(6.0*ZgradCells[celli]*deltaCells[celli]) ;
    // }
    else {
    	ZvarCells[celli] = 0.0;
    }
}

//- Updating instantaneous and effective stoichiometric dissipation rate

scalarField& chi_stCells = chi_stInstan.ref();    // instantaneous chi_st
chi_stCells=C_chi*turbulence->muEff()/(rho)*magSqr(fvc::grad(Z));

forAll(chi_st, celli){
    if (mesh.C()[celli].component(vector::Z) < 0.003){
        chi_st[celli] = chi_stCells[celli];
    }
 }

{
    //- Solving transport equation of effective chi_st
    fvScalarMatrix chi_stEqn
    (
        (
        fvm::ddt(rho, chi_st)
          + fvm::div(phi, chi_st)
          - fvm::laplacian(turbulence->muEff(), chi_st)
          -rho*2.0/0.024781*chi_st*(C_chi*turbulence->muEff()/(rho)*magSqr(fvc::grad(Z))-chi_st)
        )
    );
    chi_stEqn.relax();
    chi_stEqn.solve(mesh.solver("chi_st"));
    chi_st.correctBoundaryConditions();
}


forAll(IntermF, celli)
{
    if (Z[celli] > 0.0 && Z[celli] < 1.0)
    {
        IntermF[celli] = Zvar[celli]/(Z[celli]*(1-Z[celli]));
    }
    else 
    {
        IntermF[celli] = 0;
    }
}

