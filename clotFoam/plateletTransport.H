/**************************************************************************
    Transport the Mobile Platelets via Advection and Diffusion
***************************************************************************/
{    
    // Interpolate Theta_T to the faces for advective flux
    // Best to use interpolate(Theta_Tfa)  downwind phi;
    surfaceScalarField Theta_Tf
    (
        fvc::interpolate(Theta_T,"interpolate(Theta_Tfa)")
    );

    // Calculate the hindered transport function on the faces
    surfaceScalarField W
    ( 
        max
        (
            0*Theta_Tf,
            min
            ( 
                1.0 + 0*Theta_Tf, 
                1.0/std::tanh(std::atan(1)* 4.0 * 1.0) *
                    tanh( constant::mathematical::pi * (1.0 - Theta_Tf))
            )
        )
    );

    // Advective Flux Field
    surfaceScalarField phiPlt( phi*W );

    // Interpolate Theta_T to the faces for diffusive flux
    // Best to use interpolate(Theta_Tfa)  hinderedDiffusion or localMax;
    Theta_Tf =
        fvc::interpolate(Theta_T,"interpolate(Theta_Tfd)");

    // Calculate the hindered transport function on the faces
    W = 
        max
        (
            0*Theta_Tf,
            min
            ( 
                1.0 + 0*Theta_Tf, 
                1.0/std::tanh(std::atan(1)* 4.0 * 1.0) *
                    tanh( constant::mathematical::pi * (1.0 - Theta_Tf))
            )
        );

    // Diffusive Flux Field
    surfaceScalarField DpW( pltConst.Dp*W );

    // Transport the Mobile Platelets
    for(int i=0; i < pltConst.num_plts_mobile; i++)
    {
        fvScalarMatrix PltTransportEqn
        (
            fvm::ddt(Plt.field[i])
            + fvm::div(phiPlt,Plt.field[i],"div(phiPlt,Plt)")
            - fvm::laplacian(DpW,Plt.field[i])
        );

        PltTransportEqn.relax();
        PltTransportEqn.solve();
        Plt.field[i].correctBoundaryConditions();
    }
}

