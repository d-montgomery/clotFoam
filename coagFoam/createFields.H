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

// for openfoam-v10
//setRefCell(p, mesh.solution().dict().subDict("PISO"), pRefCell, pRefValue);
//mesh.schemes().setFluxRequired(p.name());

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

// k_adh Field
Info<< "Reading field Kadh\n" << endl;
volScalarField Kadh
(
    IOobject
    (
        "Kadh",
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
    num_plts,
    Pmax, Kadh, KcohPmax,
    eta, eta_ast, eta_t, g0, 
    ADP, kADP_Pla, ADP_ast,
    ke2_Pla, e2_ast,
    Theta_T, Theta_B
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
    Kc0_on.value(), Kc0_off.value(), Kc0_cat.value()
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
seBound.setInitialSEConc(E0_IC); // Set intial concentration on Subendothelium


/**************************************************************************
    Fluid-Phase Fields
***************************************************************************/
Info<< "Creating species object for fluidPhase Chemicals\n" << endl;
// Object for for seBound chemicals

Species_fluidPhase fluidPhase
(
    num_fluidPhase,
    N_Avo,
    Kc1_cat,
    Ke1_on, 
    Ke1_off, 
    Ks2_on, 
    Ks2_off,
    N_E1,
    N_S2
);

// Set ALL ptrList Sizes
fluidPhase.setPtrListSize();

// S1 (FX)
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

// S2
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

// E1 (prothrombinase-like)
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
    N_Avo,
    Kc1_on, Kc1_off, Kc1_cat,
    Ke1_on, Ke1_off,
    Ks2_on, Ks2_off,
    N_E1, N_S2
);

// Set ALL ptrList Sizes
pltBound.setPtrListSize();

// S2b (Prothrombin?)
pltBound.field.set
(
    0,
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

// E1b (??)
pltBound.field.set
(    
    1,
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

// C1
pltBound.field.set
(    
    2,
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

// Set all necessary fields for RK method = 0 with proper dimensions
pltBound.initializeRKFields();