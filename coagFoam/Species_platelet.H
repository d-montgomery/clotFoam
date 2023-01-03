class Species_platelet : public Species 
{    
    private:
        const dimensionedScalar& Pmax;
        const volScalarField& Kadh;
        const dimensionedScalar KcohPmax;
        const volScalarField& eta;
        const dimensionedScalar& eta_ast;
        const dimensionedScalar& eta_t;
        const dimensionedScalar& g0;
        const volScalarField& ADP;
        const dimensionedScalar& kADP_Pla;
        const dimensionedScalar& ADP_ast;
        const dimensionedScalar& ke2_Pla;
        const dimensionedScalar& e2_ast;
        volScalarField& Theta_T;
        volScalarField& Theta_B;
        Species *fluidPhase {nullptr};

    public:
        //- Constructor
        Species_platelet
        (
            const int numSpecies_,
            const dimensionedScalar& Pmax_,
            const volScalarField& Kadh_,
            const dimensionedScalar KcohPmax_,
            const volScalarField& eta_,
            const dimensionedScalar& eta_ast_,
            const dimensionedScalar& eta_t_,     
            const dimensionedScalar& g0_,
            const volScalarField& ADP_,
            const dimensionedScalar& kADP_Pla_,
            const dimensionedScalar& ADP_ast_,
            const dimensionedScalar& ke2_Pla_,
            const dimensionedScalar& e2_ast_,
            volScalarField& Theta_T_,
            volScalarField& Theta_B_
        )
        :
            Species(numSpecies_),
            Pmax(Pmax_), Kadh(Kadh_), KcohPmax(KcohPmax_),
            eta(eta_), eta_ast(eta_ast_), eta_t(eta_t_), g0(g0_), 
            ADP(ADP_), kADP_Pla(kADP_Pla_), ADP_ast(ADP_ast_),
            ke2_Pla(ke2_Pla_), e2_ast(e2_ast_),
            Theta_T(Theta_T_), Theta_B(Theta_B_)
        {}
        
        //- Member Functions
        void setPtrs(Species& fluidPhase_)
        {
            fluidPhase = &fluidPhase_;
        }
        
        // update the total and bound platelet fractions
        void updateFractions()
        {
            Theta_T = (field[0] + field[1] + field[2] + field[3])/Pmax;
            Theta_B = (field[2] + field[3])/Pmax;
        }

        // Activation hill function for any concentration C, given rate kC_Pla & C_ast
        volScalarField Kact
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
            return  g0 * max(
                            0*eta,
                            min
                            (
                                1. + 0*eta,
                                pow(eta-eta_t,3.0) /
                                ( pow(eta_ast,3.0) + pow(eta - eta_t,3.0) )
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
            Key for variable references
            Plt[0] -> Pmu,   fluidPhase[0] -> S1,   
            Plt[1] -> Pma,   fluidPhase[1] -> S2,   
            Plt[2] -> Pba,   fluidPhase[2] -> E1,   
            Plt[3] -> Pbse   fluidPhase[3] -> E2 (thrombin)

            input refers to the function input from RK4, 
            e.g. input[0] = Pmu_0 + h/2 * k1

            Output is k[i] = R_i(input)
            \* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

            const volScalarField& Pmu = input[0];
            const volScalarField& Pma = input[1];
            const volScalarField& Pba = input[2];
            const volScalarField& Pbse = input[3];
            const volScalarField& E2 = fluidPhase->field[3];
            
            // Mobile Unactivated
            k[0] = -1.*Kadh * max(0*Pbse, min(Pmax-Pbse, Pmax+0*Pbse)) * Pmu 
              -( Kact(kADP_Pla,ADP_ast,ADP) + Kact(ke2_Pla,e2_ast,E2) )* Pmu;
            
            // Mobile Activated
            k[1] = -1.*Kadh * max(0*Pbse, min(Pmax-Pbse, Pmax+0*Pbse)) * Pma 
                +( Kact(kADP_Pla,ADP_ast,ADP) + Kact(ke2_Pla,e2_ast,E2) )* Pmu
                - KcohPmax*g()*Pma;
            
            // Bound Activated
            k[2] = -1. * Kadh * max( 0*Pbse, min(Pmax-Pbse, Pmax+0*Pbse) ) * Pba
                    + KcohPmax*g()*Pma;
            
            // SE-Bound Activated
            k[3] = Kadh * max( 0*Pbse, min(Pmax-Pbse, Pmax+0*Pbse) ) 
                    * (Pmu + Pma + Pba);

        }
};
