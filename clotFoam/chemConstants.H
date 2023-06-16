struct chemConstants
{
    // Reference to the inputParameters dict
    IOdictionary& inputParameters;

    // Avogadro's Constant
    dimensionedScalar N_Avo
    {
        "N_Avo", 
        dimless/dimMoles,
        inputParameters.lookup("N_Avo")
    };

    // Number of binding sites on platelet surface for specific chems
    scalar N_1 
    {
        readScalar(inputParameters.lookup("N_1"))
    };
    scalar N_2 
    {
        readScalar(inputParameters.lookup("N_2"))
    };

    // Kinetic Rates for C0
    dimensionedScalar kC0_on
    {
        "kC0_on",
        dimVolume/dimMoles/dimTime,
        inputParameters.lookup("kC0_on")
    }; 
    dimensionedScalar kC0_off
    {
        "kC0_off",
        dimless/dimTime,
        inputParameters.lookup("kC0_off")
    }; 

    dimensionedScalar kC0_cat
    {
        "kC0_cat",
        dimless/dimTime,
        inputParameters.lookup("kC0_cat")
    }; 

    // Kinetic Rates for C1
    dimensionedScalar kC1_on
    {
        "kC1_on",
        dimVolume/dimMoles/dimTime,
        inputParameters.lookup("kC1_on")
    }; 
    dimensionedScalar kC1_off
    {
        "kC1_off",
        dimless/dimTime,
        inputParameters.lookup("kC1_off")
    }; 

    dimensionedScalar kC1_cat
    {
        "kC1_cat",
        dimless/dimTime,
        inputParameters.lookup("kC1_cat")
    }; 

    // Kinetic Rates for C2
    dimensionedScalar kC2_on
    {
        "kC2_on",
        dimVolume/dimMoles/dimTime,
        inputParameters.lookup("kC2_on")
    }; 
    dimensionedScalar kC2_off
    {
        "kC2_off",
        dimless/dimTime,
        inputParameters.lookup("kC2_off")
    }; 

    dimensionedScalar kC2_cat
    {
        "kC2_cat",
        dimless/dimTime,
        inputParameters.lookup("kC2_cat")
    }; 

    // Kinetic Rates for S1
    dimensionedScalar kS1_on
    {
        "kS1_on",
        dimVolume/dimMoles/dimTime,
        inputParameters.lookup("kS1_on")
    }; 
    dimensionedScalar kS1_off
    {
        "kS1_off",
        dimless/dimTime,
        inputParameters.lookup("kS1_off")
    }; 

    // Kinetic Rates for S2
    dimensionedScalar kS2_on
    {
        "kS2_on",
        dimVolume/dimMoles/dimTime,
        inputParameters.lookup("kS2_on")
    }; 
    dimensionedScalar kS2_off
    {
        "kS2_off",
        dimless/dimTime,
        inputParameters.lookup("kS2_off")
    }; 

    // Kinetic Rates for E1
    dimensionedScalar kE1_on
    {
        "kE1_on",
        dimVolume/dimMoles/dimTime,
        inputParameters.lookup("kE1_on")
    }; 
    dimensionedScalar kE1_off
    {
        "kE1_off",
        dimless/dimTime,
        inputParameters.lookup("kE1_off")
    }; 

    // Kinetic Rates for E2
    dimensionedScalar kE2_on
    {
        "kE2_on",
        dimVolume/dimMoles/dimTime,
        inputParameters.lookup("kE2_on")
    }; 

    dimensionedScalar kE2_off
    {
        "kE2_off",
        dimless/dimTime,
        inputParameters.lookup("kE2_off")
    }; 

    /***********************************************************************
        Constructor
    ************************************************************************/
    chemConstants(IOdictionary& inputParameters_)
    :
    inputParameters(inputParameters_)
    {}
};