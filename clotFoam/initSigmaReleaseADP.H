/******************************************************************************
    Initialize components necessary for sigma_release in ADP source
******************************************************************************/
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
    <<", sigma_N = "<< sigma_N <<","<<endl;
}
else
{
    Info<< "\n!!! Invalid sigma_dt  !!!" 
        << "\nremainder(sigma_Tf,sigma_dt) must be ~= 0" << endl;
    FatalError<< "\n!!! Invalid sigma_dt  !!!" 
              << "\nremainder(sigma_Tf,sigma_dt) must be ~= 0" 
              << exit(FatalError);
}

// Make sure that the startTime is divisible by sigma_dt
if
(
    remainder
    (
        runTime.controlDict().get<scalar>("startTime"),
        sigma_dt
    ) > 1e-12
)
{
    Info<< "\n!!! Invalid startTime !!!" 
    << "\nremainder(startTime,sigma_dt) must be ~= 0" << endl;
    FatalError<< "\n!!! Invalid startTime  !!!" 
    << "\nremainder(startTime,sigma_dt) must be ~= 0" 
            << exit(FatalError);
}

// Initialize pointer list for ddt_Pbnew
PtrList<volScalarField> ddt_Pbnew(sigma_N+1);

// Initialize IOList for tau using copy constructor of a scalarList
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
   scalarList(sigma_N+1, 0.)//initTau // copy constructor
);

// Initialize ddt_Pbnew_i as zeros, or read in if startTime > 0
// Also initializes tau if startTime > 0 and no tau present at startTime
forAll(ddt_Pbnew,i)
{
    scalar t0 = runTime.controlDict().get<scalar>("startTime");
    
    if(t0 < 1e-12) // If startTime = 0, initialize as zeros
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
                    IOobject::NO_READ,                                      
                    IOobject::NO_WRITE
                ),  
                0*pltConst.Pmax*Theta_B / runTime.deltaT()
            )   
        );  
    }
    else // If startTime > 0, must read in ddt_Pbnew[i]
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
                    IOobject::MUST_READ,                                      
                    IOobject::NO_WRITE
                ),
                mesh
            )   
        );  

        // Initialize tau
        if (tau[i] < 1e-12)
        {
            // This initializes tau properly if there is no tau read in 
            // at startTime.  This assumes that mod(startTime, sigma_dt) = 0.
            // This allows users to restart parallel codes where tau may
            // have been lost after executing reconstructPar. Some accuracy
            // may be lost if startTime > 1000. Best practice is to read in
            // tau from file, but this works if it's not present.

            tau[i] = std::max(0., t0 - sigma_Tf + sigma_dt*i );
        }
    }
 }

 // Print the first and last entries for tau
Info<<"                    tau[0] = "<< tau[0]
<< ", tau[end-1] = "<< tau[sigma_N-1] 
<< ", tau[end] = "<< tau[sigma_N] << endl; 
