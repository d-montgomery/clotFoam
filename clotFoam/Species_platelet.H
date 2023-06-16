class Species_platelet : public Species 
{    
    private:
        const plateletConstants& c;
        const volScalarField& Hadh;
        const volScalarField& eta;
        const volScalarField& ADP;
        volScalarField& Theta_T;
        volScalarField& Theta_B;
        Species *fluidPhase {nullptr};

    public:
        //- Constructor
        Species_platelet
        (
            const int numSpecies_,
            const plateletConstants& c_,
            const volScalarField& Hadh_,
            const volScalarField& eta_,
            const volScalarField& ADP_,
            volScalarField& Theta_T_,
            volScalarField& Theta_B_
        )
        :
            Species(numSpecies_),
            c(c_),
            Hadh(Hadh_), 
            eta(eta_), 
            ADP(ADP_), 
            Theta_T(Theta_T_), 
            Theta_B(Theta_B_),
            fluidPhase(nullptr)
        {}
        
        //- Member Functions
        void setPtrs(Species& fluidPhase_)
        {
            fluidPhase = &fluidPhase_;
        }
        
        // update the total and bound platelet fractions
        void updateFractions()
        {
            Theta_T = (field[0] + field[1] + field[2] + field[3])/c.Pmax;
            Theta_B = (field[2] + field[3])/c.Pmax;
        }

        // Activation hill function for any concentration C, given rate kC_Pla & C_ast
        volScalarField A
        (
            const dimensionedScalar& kC_Pla, 
            const dimensionedScalar& C_ast, 
            const volScalarField& C // (ADP or e2)
        )
        {
            return kC_Pla * max(0.*C, C) / (C_ast + max(0.*C,C) );
        }

        // Platelet binding affinity function g(eta)  
        volScalarField g()
        {
            return  c.g0 * max(
                            0*eta,
                            min
                            (
                                1. + 0*eta,
                                Foam::pow(eta - c.eta_t,3.0) /
                                ( Foam::pow(c.eta_ast,3.0) + Foam::pow(eta - c.eta_t,3.0) )
                            )
                            );
        }

        // Function that edits the k's for RK4 method.  This assumes that 
        // all k's are defined as 
        //      k[i] = R_i(input) 
        // where R represents the reactions in the ADR equation 
        //      d/dt(C_i) = - div(u C_i - D grad(C_i)) + R_i(C1,C2,...)
        void updateKs()
        {
            /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
            Key for variable references (note mobile platetlets come first)
            Plt[0] -> Pmu,   fluidPhase[0] -> S1,   
            Plt[1] -> Pma,   fluidPhase[1] -> S2,   
            Plt[2] -> Pba,   fluidPhase[2] -> E1,   
            Plt[3] -> Pbse   fluidPhase[3] -> E2 (thrombin)

            input refers to the function input from RK4, 
            e.g. input[0] = Pmu_0 + h/2 * k1

            Output is k[i] = R_i(input)
            \* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

            // Quick references to various fields
            const volScalarField& Pmu = input[0];
            const volScalarField& Pma = input[1];
            const volScalarField& Pba = input[2];
            const volScalarField& Pbse = input[3];
            const volScalarField& E2 = fluidPhase->field[3];

            // Calculate kAdh*Hadh once
            const volScalarField kAdh_x_Hadh = c.kAdh * Hadh;
            
            // Mobile Unactivated
            k[0] = -kAdh_x_Hadh * max(0*Pbse, min(c.Pmax-Pbse, c.Pmax+0*Pbse)) * Pmu 
              -( A(c.kADP_Pla,c.ADP_ast,ADP) + A(c.kE2_Pla,c.E2_ast,E2) )* Pmu;
            
            // Mobile Activated
            k[1] = -kAdh_x_Hadh * max(0*Pbse, min(c.Pmax-Pbse, c.Pmax+0*Pbse)) * Pma 
                +( A(c.kADP_Pla,c.ADP_ast,ADP) + A(c.kE2_Pla,c.E2_ast,E2) )* Pmu
                - c.kCohPmax*g()*Pma;
            
            // Bound Activated
            k[2] = -kAdh_x_Hadh * max( 0*Pbse, min(c.Pmax-Pbse, c.Pmax+0*Pbse) ) * Pba
                    + c.kCohPmax*g()*Pma;
            
            // SE-Bound Activated
            k[3] = kAdh_x_Hadh * max( 0*Pbse, min(c.Pmax-Pbse, c.Pmax+0*Pbse) ) 
                    * (Pmu + Pma + Pba);

        }
};
