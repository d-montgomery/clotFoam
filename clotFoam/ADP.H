/**************************************************************************
    Solve for ADP
***************************************************************************/
// Transport the ADP
fvScalarMatrix ADPEqn
(
    fvm::ddt(ADP)
    + fvm::div(phi,ADP,"div(phi,chems)")
    - fvm::laplacian(D_ADP,ADP)
    - sigma_release
);

ADPEqn.relax();
ADPEqn.solve();
ADP.correctBoundaryConditions();

// Update the sigma_release function
tcount = tcount + runTime.deltaTValue(); 
if ( mag(tcount - sigma_dt) <= runTime.deltaTValue() )
{
    tcount = 0.0;

    // Move previous ddt_Pbnew's and tau's down one entry
    for(int k = 0; k < sigma_N; k++)
    {
        ddt_Pbnew[k] = ddt_Pbnew[k+1]; 
        tau[k] = tau[k+1];
    }

    // Update current ddt_Pbnew and tau
    ddt_Pbnew[sigma_N] = 
        pltConst.Pmax*(Theta_B-Theta_B.oldTime()) / runTime.deltaT();
    tau[sigma_N] = runTime.time().value();

    // Integrate sigma_release
    sigma_release = 0.0*sigma_release;

    for(int k = 0; k < sigma_N; k++)
    {
        sigma_release = sigma_release + 0.5 * ADP_rel * (tau[k+1]-tau[k])
            *( R_ADP(runTime.time().value()-tau[k]) * ddt_Pbnew[k] 
            + R_ADP(runTime.time().value()-tau[k+1]) * ddt_Pbnew[k+1] );
    }   
}


