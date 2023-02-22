/**************************************************************************
    Initialize components necessary for sigma_release in ADP source
***************************************************************************/
Info << "\nInitializing the ADP source term sigma_release" << endl;

// Setup for sigma release field
scalar tcount = 0.0; // used to determine when sigma_release updates
label sigma_N = 24; // Default number of ddt_Pbnew stored 

// Make sure sigma_Tf is divisible by provided sigma_dt
if (remainder(sigma_Tf, sigma_dt) < 1e-12 && sigma_dt > 0.0)  
{
    sigma_N = round(sigma_Tf / sigma_dt);
    Info<< "\nSigma_release info: sigma_Tf = "<< sigma_Tf
    <<", sigma_dt = "<< sigma_dt
    <<", sigma_N = "<< sigma_N <<endl;
}
else
{
    Info<< "\n!!! Invalid sigma_dt  !!!" 
        << "\nremainder(sigma_Tf,sigma_dt) must be ~= 0" << endl;
    FatalError<< "\n!!! Invalid sigma_dt  !!!" 
              << "\nremainder(sigma_Tf,sigma_dt) must be ~= 0" 
              << exit(FatalError);
}

// Initialize pointer list for ddt_Pbnew and IOList for tau.
PtrList<volScalarField> ddt_Pbnew(sigma_N+1);
IOList<scalar> tau
(
    IOobject
                (   
                    "tau",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,                                      
                    IOobject::NO_WRITE
                ),  
    sigma_N+1 //num elements in List
);

// Initialize ddt_Pbnew_i as zeros, or read if present 
forAll(ddt_Pbnew,i)
{
        ddt_Pbnew.set
        (   
            i,  
            new volScalarField
            (   
                IOobject
                (   
                    "ddt_Pbnew_" + std::to_string(i),
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,                                      
                    IOobject::NO_WRITE
                ),  
                0*Pmax*Theta_B / runTime.deltaT()
            )   
        );  

        
        // Initialize tau to zero if startTime = 0
        if (runTime.controlDict().lookup<double>("startTime") < 1e-12)
        {
            tau[i] = 0.;
        }
        
}