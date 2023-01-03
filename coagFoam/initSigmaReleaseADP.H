/**************************************************************************
    Initialize components necessary for sigma_release in ADP source
***************************************************************************/
// Setup for sigma release field
scalar tcount = 0.0; // used to determine when sigma_release updates
label sigma_N = 24; // Default number of ddt_Pbnew stored 

// Make sure sigma_Tf is divisible by provided sigma_dt
if (fmod(sigma_Tf, sigma_dt) == 0.0 && sigma_dt > 0.0)  
{
    sigma_N = round(sigma_Tf / sigma_dt);
    Info<< "\nSigma_release info: sigma_Tf = "<< sigma_Tf
    <<", sigma_dt = "<< sigma_dt
    <<", sigma_N = "<< sigma_N <<endl;
}
else
{
    Info<< "\n!!! Invalid sigma_dt  !!!" <<endl;
    Info<< "fmod(sigma_Tf,sigma_dt) must be 0" << endl;
    FatalError<< nl << exit(FatalError);
}

// Initialize pointer list for ddt_Pbnew and list for tau.
PtrList<volScalarField> ddt_Pbnew(sigma_N+1);
double tau [sigma_N+1] = {};

forAll(ddt_Pbnew,i)
{
    // Initialize ddt_Pbnew as zeros, or read if present 
    // (Only reads last field!  
    //  Needs to be updated to read in all fields)
    ddt_Pbnew.set
    (   
        i,  
        new volScalarField
        (   
            IOobject
            (   
                "ddt_Pbnew",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,                                      
                IOobject::AUTO_WRITE
            ),  
            0*Pmax*Theta_B / runTime.deltaT()
        )   
    );  

    // Iinitialize tau with all zeros if startTime = 0.
    // else: tau = [ startTime - sigma_Tf, startTime - sigma_Tf+sigma_dt ,
    //               ..., startTime - sigma_dt, startTime]
    // !! This assumes startTime is an integer !!

    tau[i] = runTime.controlDict().lookupOrDefault<double>("startTime",0.0);
    /*std::max
        (
            0.0,
            runTime.controlDict().lookupOrDefault<double>("startTime",0.0)
            - sigma_Tf + sigma_dt*i
        );
        */
}