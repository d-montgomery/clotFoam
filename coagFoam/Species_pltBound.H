class Species_pltBound : public Species
{
    private:
        const dimensionedScalar& N_Avo;
        const dimensionedScalar& Kc1_on;
        const dimensionedScalar& Kc1_off; 
        const dimensionedScalar& Kc1_cat;
        const dimensionedScalar& Ke1_on;
        const dimensionedScalar& Ke1_off; 
        const dimensionedScalar& Ks2_on;
        const dimensionedScalar& Ks2_off; 
        const scalar& N_E1; 
        const scalar& N_S2; 
        Species *Plt {nullptr};
        Species *fluidPhase {nullptr};
        

    public: 
        //- Constructor
        Species_pltBound
        (
            const int& numSpecies_,
            const dimensionedScalar& N_Avo_,
            const dimensionedScalar& Kc1_on_, 
            const dimensionedScalar& Kc1_off_, 
            const dimensionedScalar& Kc1_cat_,
            const dimensionedScalar& Ke1_on_, 
            const dimensionedScalar& Ke1_off_,
            const dimensionedScalar& Ks2_on_, 
            const dimensionedScalar& Ks2_off_,
            const scalar& N_E1_, 
            const scalar& N_S2_
        ) 
        : 
            Species(numSpecies_), 
            N_Avo(N_Avo_),
            Kc1_on(Kc1_on_), Kc1_off(Kc1_off_), Kc1_cat(Kc1_cat_),
            Ke1_on(Ke1_on_), Ke1_off(Ke1_off_),
            Ks2_on(Ks2_on_), Ks2_off(Ks2_off_),
            N_E1(N_E1_), N_S2(N_S2_)
        {}

        void setPtrs(Species& Plt_, Species& fluidPhase_)
        {
            Plt = &Plt_;
            fluidPhase = &fluidPhase_;
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
            pltBound[0] = S2b,     fluidPhase[0] -> S1,     Plt[0] -> Pmu,
            pltBound[1] = E1b,     fluidPhase[1] -> S2,     Plt[1] -> Pma,
            pltBound[2] = C1,      fluidPhase[2] -> E1,     Plt[2] -> Pba,
                                   fluidPhase[3] -> E2,     Plt[3] -> Pbse.
                                    
            input refers to the function input from RK4, 
            e.g. input[0] = S2_b + h/2 * k1

            Output is k[i] = R_i(input)
            \* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

            // pltBound Reference Vars
            const volScalarField& S2b = input[0];
            const volScalarField& E1b = input[1];
            const volScalarField& C1 = input[2];

             // fluidPhase Reference Vars
            //const volScalarField& S1 = fluidPhase->input[0];
            const volScalarField& S2 = fluidPhase->input[1];
            const volScalarField& E1 = fluidPhase->input[2];

            // Platelet Reference Vars
            const volScalarField& Pba = Plt->input[2];
            const volScalarField& Pbse = Plt->input[3];

            // S2b
            k[0] = Ks2_on * N_S2 * (Pba + Pbse) * S2 / N_Avo
                   -Ks2_off * S2b 
                   -Kc1_on * S2b * E1b 
                   +Kc1_off*C1;
            
            // E1b
            k[1] = Ke1_on * N_E1 * (Pba + Pbse) * E1 / N_Avo
                   -Ke1_off * E1b
                   -Kc1_on * S2b * E1b
                   +(Kc1_off + Kc1_cat)*C1;
            
            // C1
            k[2] = Kc1_on * S2b * E1b 
                   -(Kc1_off + Kc1_cat)*C1;
        }
};