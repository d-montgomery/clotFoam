/**************************************************************************
    Fluid Fields
***************************************************************************/
// Pressure Field
Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Velocity Field
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

/**************************************************************************
    Fields for Platelet Reactions
***************************************************************************/
// Virtual substance secreted from Bound Platelets
Info<< "Reading field eta\n" << endl;
volScalarField eta
(
    IOobject
    (
        "eta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Hadh Field (defines adhesion region)
Info<< "Reading field Hadh\n" << endl;
volScalarField Hadh
(
    IOobject
    (
        "Hadh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

/**************************************************************************
    Fields for ADP and ADP source term
***************************************************************************/
// ADP Field
Info<< "Reading field ADP\n" << endl;
volScalarField ADP
(
    IOobject
    (
        "ADP",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// sigma_release Field
Info<< "Reading field sigma_release\n" << endl;
volScalarField sigma_release
(
    IOobject
    (
        "sigma_release",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

/**************************************************************************
    Platelet Fraction Fields
***************************************************************************/
// Total Platelet Fraction
Info<< "Reading field Theta_T\n" << endl;
volScalarField Theta_T
(
    IOobject
    (
        "Theta_T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Bound Platelet Fraction
Info<< "Reading field Theta_B\n" << endl;
volScalarField Theta_B
(
    IOobject
    (
        "Theta_B",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

/**************************************************************************
    Platelet Fields
***************************************************************************/

// Object for for platelets
Species_platelet Plt
(
    pltConst.num_plts,
    pltConst,
    Hadh, 
    eta,  
    ADP, 
    Theta_T, 
    Theta_B
);

// Set ALL ptrList Sizes
Plt.setPtrListSize();

// Plt.field[0] -> Mobile Unactivated
// Plt.field[1] -> Mobile Activated
// Plt.field[2] -> Bound Activated
// Plt.field[3] -> SEBound Activated

// Mobile Unactivated Platelet Field
Plt.field.set
(
    0,
    new volScalarField
    (   
        IOobject
        (   
            "plt_Pmu",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,                                      
            IOobject::AUTO_WRITE
        ),  
        mesh
    )   
);

// Mobile Activated Platelet Field
Plt.field.set
(    
    1,
    new volScalarField
    (
        IOobject
        (
            "plt_Pma",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);

// Bound Activated Platelet Field
Info<< "Reading field Pba\n" << endl;
Plt.field.set
(
    2,
    new volScalarField
    (
        IOobject
        (
            "plt_Pba",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);

// Subendothelial-Bound Platelet Field
Info<< "Reading field Pbse\n" << endl;
Plt.field.set
(
    3, 
    new volScalarField
    (
        IOobject
        (
            "plt_Pbse",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);

// Set all necessary fields for RK method = 0*Plt.field
Plt.initializeRKFields();

/**************************************************************************
    SE-Bound Fields
***************************************************************************/
Info<< "Creating species object for seBound Chemicals\n" << endl;

// Object for for seBound chemicals
Species_seBound seBound
(
    num_seBound,
    injuryPatchID,
    chemConst
);

// Set ALL ptrList Sizes
seBound.setPtrListSize();

// E0 (TF:VIIa)
seBound.field.set
(
    0,
    new volScalarField
    (   
        IOobject
        (   
            "seBound_E0",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,                                      
            IOobject::AUTO_WRITE
        ),  
        mesh
    )   
);

// C0 
seBound.field.set
(    
    1,
    new volScalarField
    (
        IOobject
        (
            "seBound_C0",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);

// Set all necessary fields for RK method = 0 with proper dimensions
seBound.initializeRKFields();

// Set the initial concentration for the seBound species
seBound.setInitialConcentration();

/**************************************************************************
    Fluid-Phase Fields
***************************************************************************/
Info<< "Creating species object for fluidPhase Chemicals\n" << endl;
// Object for for seBound chemicals

Species_fluidPhase fluidPhase
(
    num_fluidPhase,
    chemConst
);

// Set ALL ptrList Sizes
fluidPhase.setPtrListSize();

// S1 (FV and FX)
fluidPhase.field.set
(
    0,
    new volScalarField
    (   
        IOobject
        (   
            "fluidPhase_S1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,                                      
            IOobject::AUTO_WRITE
        ),  
        mesh
    )   
);

// S2 (Prothrombin)
fluidPhase.field.set
(
    1,
    new volScalarField
    (   
        IOobject
        (   
            "fluidPhase_S2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,                                      
            IOobject::AUTO_WRITE
        ),  
        mesh
    )   
);

// E1 (FXa)
fluidPhase.field.set
(    
    2,
    new volScalarField
    (
        IOobject
        (
            "fluidPhase_E1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);

// E2 (thrombin)
fluidPhase.field.set
(    
    3,
    new volScalarField
    (
        IOobject
        (
            "fluidPhase_E2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);

// Set all necessary fields for RK method = 0 with proper dimensions
fluidPhase.initializeRKFields();

/**************************************************************************
    Platelet-Bound Fields
***************************************************************************/

Info<< "Reading PtrList for pltBound\n" << endl;

Species_pltBound pltBound
(
    num_pltBound,
    chemConst
);

// Set ALL ptrList Sizes
pltBound.setPtrListSize();

// S1b (FV)
pltBound.field.set
(
    0,
    new volScalarField
    (   
        IOobject
        (   
            "pltBound_S1b",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,                                      
            IOobject::AUTO_WRITE
        ),  
        mesh
    )   
);

// S2b (Prothrombin)
pltBound.field.set
(
    1,
    new volScalarField
    (   
        IOobject
        (   
            "pltBound_S2b",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,                                      
            IOobject::AUTO_WRITE
        ),  
        mesh
    )   
);

// E1b (Prothrombinase)
pltBound.field.set
(    
    2,
    new volScalarField
    (
        IOobject
        (
            "pltBound_E1b",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);

// E2b (Bound Thrombin)
pltBound.field.set
(    
    3,
    new volScalarField
    (
        IOobject
        (
            "pltBound_E2b",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);

// C1 (S2b:E1b)
pltBound.field.set
(    
    4,
    new volScalarField
    (
        IOobject
        (
            "pltBound_C1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);

// C2 (S1b:E2b)
pltBound.field.set
(    
    5,
    new volScalarField
    (
        IOobject
        (
            "pltBound_C2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);

// Set all necessary fields for RK method = 0 with proper dimensions
pltBound.initializeRKFields();