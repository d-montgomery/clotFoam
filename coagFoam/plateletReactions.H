/**************************************************************************
    React the Platelets M times
***************************************************************************/
{    
    // Solve for substance eta
    eta = max(0*Theta_B, min(1.+0*Theta_B,Theta_B));
    eta.correctBoundaryConditions(); 
    solve(
        fvm::Sp(1.0,eta)
        - fvc::Sp(1.0,eta)
        -fvm::laplacian(l_D_eta*l_D_eta/4.0,eta)
    ); 
    eta.correctBoundaryConditions(); 
    eta = max(0*eta, min(1.+0*eta,eta));

    // Solve for platelet reactions M_rxn times
    for (int m = 0; m < M_rxn; m++ )
    {      
        //- RK4 Solver
        odeSolver(h_rxn, ptrToPlt, num_pltObjs);
    }
}

// Update the BC's for the mobile Platelet Species
for(int i=0; i < num_plts_mobile; i++)
{
    Plt.field[i].correctBoundaryConditions();
}
