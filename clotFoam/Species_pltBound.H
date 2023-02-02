class Species_pltBound : public Species
{
    private:
        const dimensionedScalar& N_Avo;
        const dimensionedScalar& Kc1_on;
        const dimensionedScalar& Kc1_off; 
        const dimensionedScalar& Kc1_cat;
        const dimensionedScalar& Kc2_on;
        const dimensionedScalar& Kc2_off; 
        const dimensionedScalar& Kc2_cat;
        const dimensionedScalar& Ke1_on;
        const dimensionedScalar& Ke1_off; 
        const dimensionedScalar& Ke2_on;
        const dimensionedScalar& Ke2_off; 
        const dimensionedScalar& Ks1_on;
        const dimensionedScalar& Ks1_off;
        const dimensionedScalar& Ks2_on;
        const dimensionedScalar& Ks2_off; 
        const scalar& N_1; 
        const scalar& N_2; 
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
            const dimensionedScalar& Kc2_on_, 
            const dimensionedScalar& Kc2_off_, 
            const dimensionedScalar& Kc2_cat_,
            const dimensionedScalar& Ke1_on_, 
            const dimensionedScalar& Ke1_off_,
            const dimensionedScalar& Ke2_on_, 
            const dimensionedScalar& Ke2_off_,
            const dimensionedScalar& Ks1_on_, 
            const dimensionedScalar& Ks1_off_,
            const dimensionedScalar& Ks2_on_, 
            const dimensionedScalar& Ks2_off_,
            const scalar& N_1_, 
            const scalar& N_2_
        ) 
        : 
            Species(numSpecies_), 
            N_Avo(N_Avo_),
            Kc1_on(Kc1_on_), Kc1_off(Kc1_off_), Kc1_cat(Kc1_cat_),
            Kc2_on(Kc2_on_), Kc2_off(Kc2_off_), Kc2_cat(Kc2_cat_),
            Ke1_on(Ke1_on_), Ke1_off(Ke1_off_),
            Ke2_on(Ke2_on_), Ke2_off(Ke2_off_),
            Ks1_on(Ks1_on_), Ks1_off(Ks1_off_),
            Ks2_on(Ks2_on_), Ks2_off(Ks2_off_),
            N_1(N_1_), N_2(N_2_)
        {}

        void setPtrs(Species& Plt_, Species& fluidPhase_)
        {
            Plt = &Plt_;
            fluidPhase = &fluidPhase_;
        }

        // Function that edits the k's for RK4 method.
        void updateKs()
        {
            /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
            This assumes that all k's are defined as 
                k[i] = R_i(input) 
            where R represents the reactions in the ADR equation 
                d/dt(C_i) = - div(u C_i - D grad(C_i)) + R_i(C1,C2,...)
            
            Key for variable references:
            pltBound[0] = S1b,     fluidPhase[0] -> S1,     Plt[0] -> Pmu,
            pltBound[1] = S2b,     fluidPhase[1] -> S2,     Plt[1] -> Pma,
            pltBound[2] = E1b,     fluidPhase[2] -> E1,     Plt[2] -> Pba,
            pltBound[3] = E2b,     fluidPhase[3] -> E2,     Plt[3] -> Pbse.
            pltBound[4] = C1,
            pltBound[5] = C2,
            
            input refers to the function input from RK4, 
            e.g. input[0] = S2_b + h/2 * k1

            Output is k[i] = R_i(input)
            \* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

            // pltBound Reference Vars
            const volScalarField& S1b = input[0];
            const volScalarField& S2b = input[1];
            const volScalarField& E1b = input[2];
            const volScalarField& E2b = input[3];
            const volScalarField& C1 = input[4];
            const volScalarField& C2 = input[5];

             // fluidPhase Reference Vars
            const volScalarField& S1 = fluidPhase->input[0];
            const volScalarField& S2 = fluidPhase->input[1];
            const volScalarField& E1 = fluidPhase->input[2];
            const volScalarField& E2 = fluidPhase->input[3];

            // Platelet Reference Vars
            const volScalarField& Pba = Plt->input[2];
            const volScalarField& Pbse = Plt->input[3];

            // S1b
            k[0] = Ks1_on * (N_1 * (Pba + Pbse) / N_Avo
                   - (S1b + E1b + C1 + C2)) * S1
                   - Ks1_off * S1b
                   - Kc2_on * S1b * E2b
                   + Kc2_off * C2;

            // S2b
            k[1] = Ks2_on * (N_2 * (Pba + Pbse) / N_Avo 
                   - (S2b + E2b + C1 + C2))* S2 
                   - Ks2_off * S2b 
                   - Kc1_on * S2b * E1b 
                   + Kc1_off*C1;
            
            // E1b
            k[2] = Ke1_on * (N_1 * (Pba + Pbse) / N_Avo 
                   - (S1b + E1b + C1 + C2))* E1 
                   - Ke1_off * E1b
                   - Kc1_on * S2b * E1b
                   + (Kc1_off + Kc1_cat)*C1
                   + Kc2_cat*C2;
            
            // E2b
            k[3] = Ke2_on * (N_2 * (Pba + Pbse) / N_Avo 
                   - (S2b + E2b + C1 + C2))* E2
                   - Ke2_off * E2b
                   - Kc2_on * S1b * E2b
                   + (Kc2_off + Kc2_cat)*C2
                   + Kc1_cat * C1;

            // C1
            k[4] = Kc1_on * S2b * E1b 
                   -(Kc1_off + Kc1_cat)*C1;

            // C2
            k[5] = Kc2_on * S1b * E2b 
                   -(Kc2_off + Kc2_cat)*C2;
        }
};