/**************************************************************************
    Transport the fluid-phase chemicals via Advection and Diffusion
***************************************************************************/
{    
    // Transport the Mobile Platelets
    forAll(fluidPhase.field, i)
    {
        fvScalarMatrix chemTransportEqn
        (
            fvm::ddt(fluidPhase.field[i])
            + fvm::div(phi,fluidPhase.field[i],"div(phi,chems)")
            - fvm::laplacian(Dc,fluidPhase.field[i])
        );

        chemTransportEqn.relax();
        chemTransportEqn.solve();
        fluidPhase.field[i].correctBoundaryConditions();
    }
}

