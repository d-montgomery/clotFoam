// Calculate the frictional resistance (Darcy Term) for NSB solver
volScalarField alpha
(
    nu*C_CK*sqr(0.6*Theta_B) / (pow(1.-0.6*Theta_B,3.))
);

// Momentum predictor
fvVectorMatrix UEqn
(
    fvm::ddt(U)
  + fvm::div(phi,U)
  - fvm::laplacian(nu,U)
  + fvm::Sp(alpha,U)
);

if (piso.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));
}

// --- PISO loop
while (piso.correct())
{
    volScalarField rAU(1.0/UEqn.A());
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
    );

    adjustPhi(phiHbyA, U, p);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, U, phiHbyA, rAU);

    // Non-orthogonal pressure corrector loop
    while (piso.correctNonOrthogonal())
    {
        // Pressure corrector

        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (piso.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

    if (runTime.write())
    {
        #include "continuityErrs.H"
    }
    U = HbyA - rAU*fvc::grad(p);
    U.correctBoundaryConditions();
}
