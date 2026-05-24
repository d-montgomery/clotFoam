Info<< "\nReading inputParameters\n" << endl;

IOdictionary inputParameters
(
    IOobject
    (
        "inputParameters",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Patch identifier of the injury boundary
label injuryPatchID = mesh.boundaryMesh().findPatchID("injuryWalls");

// Flag for whether coagulation reactions occur or not 
const bool coagReactionsOn =
            runTime.controlDict().lookupOrDefault("coagReactionsOn", true);

// Flag for whether or not to smooth the adhesion region
const bool smoothHadh =
            runTime.controlDict().lookupOrDefault("smoothHadh",false);
            
/**************************************************************************
    Fluid Scalars and Parameters
***************************************************************************/

// Dynamic Viscosity
dimensionedScalar nu
(
    "nu",
    dimViscosity,
    inputParameters.get<scalar>("nu")
);

// Carman-Kozeny Consant
dimensionedScalar C_CK
(
    "C_CK",
    dimless/dimArea,
    inputParameters.get<scalar>("C_CK")
);

/**************************************************************************
    Constants for Platelets
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 See plateletConstants.H. The constants for platelets are created in a 
 structure so they can be easily passed to various classes.  This makes 
 the code easier to modify for more complex platelet/coagulation models.
 ***************************************************************************/

// Create an instance of plateletConstants for passing to various classes
Info<< "\nCreating an instance of plateletConstants\n" << endl;
plateletConstants pltConst(inputParameters);

/**************************************************************************
    Constants for ADP
***************************************************************************/
Info<< "\nReading constants for ADP and sigma_release\n" << endl;
// ADP Diffusion Constant
dimensionedScalar D_ADP
(
    "D_ADP",
    dimViscosity,
    inputParameters.get<scalar>("D_ADP")
);

// Total ADP Released (\hat A)
dimensionedScalar ADP_rel
(
    "ADP_rel",
    dimMoles,
    inputParameters.get<scalar>("ADP_rel")
);

// How often the sigma_release function is calculated
scalar sigma_dt
(
    readScalar(inputParameters.lookup("sigma_dt"))
);

// The length of time sigma_release is calculated (0 -> sigma_Tf)
scalar sigma_Tf
(
    readScalar(inputParameters.lookup("sigma_Tf"))
);

/**************************************************************************
    Params for Chemical Reactions (not Kinetic Rates)
***************************************************************************/

// Number of sebound species
label num_seBound
(
    readLabel(inputParameters.lookup("num_seBound"))
);

// Number of sebound species
label num_fluidPhase
(
    readLabel(inputParameters.lookup("num_fluidPhase"))
);

// Number of platelet-bound chem species
label num_pltBound
(
    readLabel(inputParameters.lookup("num_pltBound"))
);

// Number of time steps for reactions, so that h_rxn = dt / M_rxn
label M_rxn
(   
    readLabel(inputParameters.lookup("M_rxn"))
);

// Time step for reaction solver
dimensionedScalar h_rxn //time step for reactions
(   "h_rxn", 
    dimTime,
    1.0
);

// Diffusion Constant for all chemicals
dimensionedScalar Dc
(
    "Dc",
    dimViscosity,
    inputParameters.get<scalar>("Dc")
);


/**************************************************************************
    Kinetic Rates and Parameters for Chemical Reactions
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 See chemConstants.H. The constants for chemical reactions are created in a 
 structure so they can be easily passed to various classes.  This makes 
 the code easier to modify for more complex platelet/coagulation models. 
***************************************************************************/
// Create an instance of chemConstants for passing to various classes
Info<< "\nCreating an instance of chemConstants\n" << endl;
chemConstants chemConst(inputParameters);

