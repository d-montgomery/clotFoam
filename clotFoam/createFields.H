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

// Shear rate field 
Info<< "Reading field shearRate\n" << endl;
volScalarField shearRate
(
    IOobject
    (
        "shearRate",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    Foam::sqrt(2.0) * mag(symm( fvc::grad(U) ))
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
if(smoothHadh)
{
    Info<< "Smoothing the Hadh field \n" << endl;
    dimensionedScalar D_Hadh = 0.25 * pltConst.Pdiam;
    Hadh = Hadh + fvc::laplacian(D_Hadh*D_Hadh/4., Hadh);
    Hadh = max(0*Hadh, min(1.+0*Hadh,Hadh));    
}

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
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    0.0*eta
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
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    0.0*eta
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
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimMoles/dimVolume, 0.),
    eta.boundaryField().types()
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
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimMoles/dimVolume/dimTime, 0.),
    eta.boundaryField().types()
);

/**************************************************************************
    Platelet Object
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
1) Create derived Species object Plt
2) Set the ptrList sizes
3) Read in IC/BC for the following fields:
        Plt.field[0] -> Mobile Unactivated
        Plt.field[1] -> Mobile Activated
        Plt.field[2] -> Bound Activated
        Plt.field[3] -> SE Bound
4) Initialize the PtrList/fields needed for RK solver
***************************************************************************/
Info<< "Creating species object for seBound Chemicals\n" << endl;

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

// Mobile Unactivated Platelet Field
Info<< "Reading field plt_Pmu\n" << endl;
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
Info<< "Reading field plt_Pma\n" << endl;
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
Info<< "Reading/creating field plt_Pba\n" << endl;
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
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        0.0*Plt.field[1]
    )
);

// Subendothelial-Bound Platelet Field
Info<< "Reading/creating field plt_Pbse\n" << endl;
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
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        0.0*Plt.field[1]
    )
);

// Set all necessary fields for RK method = 0*Plt.field
Plt.initializeRKFields();

/**************************************************************************
    SE-Bound Object
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
1) Create derived Species object seBound
2) Set the ptrList sizes
3) Read in IC/BC for the following fields:
        seBound.field[0] -> E0
        seBound.field[1] -> C0
4) Initialize the PtrList/fields needed for RK solver
5) Set initial conditions for seBound species
***************************************************************************/
Info<< "Creating species object for seBound Chemicals\n" << endl;

// Object for for seBound chemicals
Species_seBound seBound
(
    num_seBound,
    injuryPatchID,
    chemConst
);

// Initialize all seBound fields if coagReactionsOn == true
if (coagReactionsOn) 
{
    // Set ALL ptrList Sizes
    seBound.setPtrListSize();

    // E0 (TF:VIIa)
    Info<< "Reading field seBound_E0\n" << endl;
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
    Info<< "Reading field seBound_C0\n" << endl;
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
}
else 
{
    Info<< "Destroying species object for seBound Chemicals\n" << endl;
    seBound.~Species();
}

/**************************************************************************
    Fluid-Phase Object
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
1) Create derived Species object fluidPhase
2) Set the ptrList sizes
3) Read in IC/BC for the following fields:
        fluidPhase.field[0] -> S1
        fluidPhase.field[1] -> S2
        fluidPhase.field[2] -> E1
        fluidPhase.field[3] -> E2
4) Initialize the PtrList/fields needed for RK solver
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
Info<< "Reading field fluidPhase_S1\n" << endl;
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
Info<< "Reading field fluidPhase_S2\n" << endl;
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
Info<< "Reading field fluidPhase_E1\n" << endl;
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
Info<< "Reading field fluidPhase_E2\n" << endl;
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
    Platelet-Bound Object
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
1) Create derived Species object fluidPhase
2) Set the ptrList sizes
3) Read in IC/BC for the following fields:
        pltBound.field[0] -> S1b
        pltBound.field[1] -> S2b
        pltBound.field[2] -> E1b
        pltBound.field[3] -> E2b
        pltBound.field[4] -> C1
        pltBound.field[5] -> C2
4) Initialize the PtrList/fields needed for RK solver
***************************************************************************/

Info<< "Creating species object for pltBound Chemicals\n" << endl;
Species_pltBound pltBound
(
    num_pltBound,
    chemConst
);

if(coagReactionsOn)
{
    // Set ALL ptrList Sizes
    pltBound.setPtrListSize();

    // S1b (FV)
    Info<< "Reading field pltBound_S1b\n" << endl;
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
    Info<< "Reading/creating field pltBound_S2b\n" << endl;
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
                IOobject::READ_IF_PRESENT,                                      
                IOobject::AUTO_WRITE
            ),  
            0.0*pltBound.field[0]
        )   
    );

    // E1b (Prothrombinase)
    Info<< "Reading/creating field pltBound_E1b\n" << endl;
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
                IOobject::READ_IF_PRESENT,                                      
                IOobject::AUTO_WRITE
            ),  
            0.0*pltBound.field[0]
        )
    );

    // E2b (Bound Thrombin)
    Info<< "Reading/creating field pltBound_E2b\n" << endl;
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
                IOobject::READ_IF_PRESENT,                                      
                IOobject::AUTO_WRITE
            ),  
            0.0*pltBound.field[0]
        )
    );

    // C1 (S2b:E1b)
    Info<< "Reading/creating field pltBound_C1\n" << endl;
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
                IOobject::READ_IF_PRESENT,                                      
                IOobject::AUTO_WRITE
            ),  
            0.0*pltBound.field[0]
        )
    );

    // C2 (S1b:E2b)
    Info<< "Reading/creating field pltBound_C2\n" << endl;
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
                IOobject::READ_IF_PRESENT,                                      
                IOobject::AUTO_WRITE
            ),  
            0.0*pltBound.field[0]
        )
    );

    // Set all necessary fields for RK method = 0 with proper dimensions
    pltBound.initializeRKFields();
}
else 
{
    Info<< "Destroying species object for pltBound Chemicals\n" << endl;
    pltBound.~Species();
}
