Info<< "Reading inputParameters\n" << endl;

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

/**************************************************************************
    Fluid Scalars and Parameters
***************************************************************************/

// Dynamic Viscosity
dimensionedScalar nu
(
    "nu",
    dimViscosity,
    inputParameters.lookup("nu")
);

// Carman-Kozeny Consant
dimensionedScalar C_CK
(
    "C_CK",
    dimless/dimArea,
    inputParameters.lookup("C_CK")
);


/**************************************************************************
    Constants for Platelets
***************************************************************************/
// Total Number of platelet species
label num_plts
(
    readLabel(inputParameters.lookup("num_plts"))
);

// Total Number of platelet species
label num_plts_mobile
(
    readLabel(inputParameters.lookup("num_plts_mobile"))
);

// Diameter of a Platelet 
dimensionedScalar Pdiam
(
    "Pdiam",
    dimLength,
    inputParameters.lookup("Pdiam")
);

// Maximum Packing Limit of Platelets
dimensionedScalar Pmax
(
    "Pmax",
    dimless/dimVolume,
    inputParameters.lookup("Pmax")
);

// Platelet Diffusion Constant
dimensionedScalar Dp
(
    "Dp",
    dimViscosity,
    inputParameters.lookup("Dp")
);

// Platelet Cohesion binding rate
dimensionedScalar KcohPmax
(
    "KcohPmax",
    dimless/dimTime,
    inputParameters.lookup("KcohPmax")
);

// Platelet activation rate by ADP k_{adp}^{pla}
dimensionedScalar kADP_Pla
(
    "kADP_Pla",
    dimless/dimTime,
    inputParameters.lookup("kADP_Pla")
);

// Platelet activation rate by ADP k_{adp}^{pla}
dimensionedScalar ke2_Pla
(
    "ke2_Pla",
    dimless/dimTime,
    inputParameters.lookup("ke2_Pla")
);

// Critcal concentration of thrombin for platelet activation
dimensionedScalar e2_ast
(
    "e2_ast",
    dimMoles/dimVolume,
    inputParameters.lookup("e2_ast")
);

/**************************************************************************
    Constants for Virtual Substance eta
***************************************************************************/
// Virtual substance Eta threshold
dimensionedScalar eta_t
(
    "eta_t",
    dimless,
    inputParameters.lookup("eta_t")
);

// Virtual substance Eta_ast
dimensionedScalar eta_ast
(
    "eta_ast",
    dimless,
    inputParameters.lookup("eta_ast")
);


// Calculate g0 for the binding affinity function g(eta)
dimensionedScalar g0
(
    "g0",
    dimless,
    0.0
);
g0 = ( pow(eta_ast,3.0) + pow(1.0 - eta_t,3.0) ) / pow(1.0-eta_t, 3.0);

/**************************************************************************
    Constants for ADP
***************************************************************************/
// ADP Diffusion Constant
dimensionedScalar D_ADP
(
    "D_ADP",
    dimViscosity,
    inputParameters.lookup("D_ADP")
);

// Critical concentration of ADP [ADP]^\ast
dimensionedScalar ADP_ast
(
    "ADP_ast",
    dimMoles/dimVol,
    inputParameters.lookup("ADP_ast")
);

// Total ADP Released (\hat A)
dimensionedScalar ADP_rel
(
    "ADP_rel",
    dimMoles,
    inputParameters.lookup("ADP_rel")
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
// Avogadro's Constant
dimensionedScalar N_Avo
(
    "N_Avo", 
    dimless/dimMoles,
    inputParameters.lookup("N_Avo")
);

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
    inputParameters.lookup("Dc")
);

// Initial concentration of E0 on the subendothelium
scalar E0_IC = readScalar(inputParameters.lookup("E0_IC"));

// Number of binding sites on platelet surface for specific chems
scalar N_E1 = readScalar(inputParameters.lookup("N_E1"));
scalar N_S2 = readScalar(inputParameters.lookup("N_S2"));

/**************************************************************************
    Kinetic Rates for Chemical Reactions 
***************************************************************************/
// Kinetic Rates for C0
dimensionedScalar Kc0_on
(
    "Kc0_on",
    dimVolume/dimMoles/dimTime,
    inputParameters.lookup("Kc0_on")
); 
dimensionedScalar Kc0_off
(
    "Kc0_off",
    dimless/dimTime,
    inputParameters.lookup("Kc0_off")
); 

dimensionedScalar Kc0_cat
(
    "Kc0_cat",
    dimless/dimTime,
    inputParameters.lookup("Kc0_cat")
); 

// Kinetic Rates for C1
dimensionedScalar Kc1_on
(
    "Kc1_on",
    dimVolume/dimMoles/dimTime,
    inputParameters.lookup("Kc1_on")
); 
dimensionedScalar Kc1_off
(
    "Kc1_off",
    dimless/dimTime,
    inputParameters.lookup("Kc1_off")
); 

dimensionedScalar Kc1_cat
(
    "Kc1_cat",
    dimless/dimTime,
    inputParameters.lookup("Kc1_cat")
); 

// Kinetic Rates for S2
dimensionedScalar Ks2_on
(
    "Ks2_on",
    dimVolume/dimMoles/dimTime,
    inputParameters.lookup("Ks2_on")
); 
dimensionedScalar Ks2_off
(
    "Ks2_off",
    dimless/dimTime,
    inputParameters.lookup("Ks2_off")
); 

// Kinetic Rates for E1
dimensionedScalar Ke1_on
(
    "Ke1_on",
    dimVolume/dimMoles/dimTime,
    inputParameters.lookup("Ke1_on")
); 
dimensionedScalar Ke1_off
(
    "Ke1_off",
    dimless/dimTime,
    inputParameters.lookup("Ke1_off")
); 