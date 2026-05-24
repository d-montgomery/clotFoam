struct plateletConstants
{
    // Reference to the inputParameters dict
    IOdictionary& inputParameters;

    // Total number of platelets and number of mobile platelets
    label num_plts
    {
        readLabel(inputParameters.lookup("num_plts"))
    };
    label num_plts_mobile
    {
        readLabel(inputParameters.lookup("num_plts_mobile"))
    };

    // Platelet diameter length
    dimensionedScalar Pdiam
    {
        "Pdiam",
        dimLength,
        inputParameters.get<scalar>("Pdiam")
    };

    // Maximum Packing Limit of Platelets
    dimensionedScalar Pmax
    {
        "Pmax",
        dimless/dimVolume,
        inputParameters.get<scalar>("Pmax")
    };

    // Platelet Diffusion Constant
    dimensionedScalar Dp
    {
        "Dp",
        dimViscosity,
        inputParameters.get<scalar>("Dp")
    };

    // Platelet adhesion binding rate
    dimensionedScalar kAdh
    {
        "kAdh",
        dimVolume/dimTime,
        inputParameters.get<scalar>("kAdh")
    };

    // Platelet Cohesion binding rate
    dimensionedScalar kCohPmax
    {
        "kCohPmax",
        dimless/dimTime,
        inputParameters.get<scalar>("kCohPmax")
    };

    // Platelet activation rate by ADP k_{adp}^{pla}
    dimensionedScalar kADP_Pla
    {
        "kADP_Pla",
        dimless/dimTime,
        inputParameters.get<scalar>("kADP_Pla")
    };

    // Critical concentration of ADP [ADP]^\ast for platelet activation
    dimensionedScalar ADP_ast
    {
        "ADP_ast",
        dimMoles/dimVol,
        inputParameters.get<scalar>("ADP_ast")
    };

    // Platelet activation rate by thrombin k_{e2}^{pla}
    dimensionedScalar kE2_Pla
    {
        "kE2_Pla",
        dimless/dimTime,
        inputParameters.get<scalar>("kE2_Pla")
    };

    // Critcal concentration of thrombin for platelet activation
    dimensionedScalar E2_ast
    {
        "E2_ast",
        dimMoles/dimVolume,
        inputParameters.get<scalar>("E2_ast")
    };

    // Critical concentration of virtual substance eta, eta_ast
    scalar eta_ast
    {
        readScalar(inputParameters.lookup("eta_ast"))
    };

    // Threshold for virtual substance eta, eta_t
    scalar eta_t
    {
        readScalar(inputParameters.lookup("eta_t"))
    };

    // Lenght of diffusion for virtual substance eta, l_D_eta
    dimensionedScalar l_D_eta
    {
        "l_D_eta",
        dimLength,
        inputParameters.get<scalar>("l_D_eta")
    };

    // Normalization constant for binding affinity function g(eta), g0
    scalar g0
    {
        ( Foam::pow(eta_ast,3.0) 
                + Foam::pow(1.0 - eta_t,3.0) ) 
               / Foam::pow(1.0-eta_t, 3.0)
    };

    /***********************************************************************
        Constructor
    ************************************************************************/
    plateletConstants(IOdictionary& inputParameters_)
    :
    inputParameters(inputParameters_)
    {}
};
