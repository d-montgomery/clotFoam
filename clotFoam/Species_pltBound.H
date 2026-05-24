class Species_pltBound : public Species
{
    private:
        const chemConstants& c;
        Species *Plt {nullptr};
        Species *fluidPhase {nullptr};
  
    public: 
        //- Constructor
        Species_pltBound
        (
            const int& numSpecies_,
            const chemConstants& c_
        ) 
        : 
            Species(numSpecies_),
            c(c_),
            Plt(nullptr),
            fluidPhase(nullptr)
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
            k[0] = c.kS1_on * (c.N_1 * (Pba + Pbse) / c.N_Avo
                   - (S1b + E1b + C1 + C2)) * S1
                   - c.kS1_off * S1b
                   - c.kC2_on * S1b * E2b
                   + c.kC2_off * C2;

            // S2b
            k[1] = c.kS2_on * (c.N_2 * (Pba + Pbse) / c.N_Avo 
                   - (S2b + E2b + C1 + C2))* S2 
                   - c.kS2_off * S2b 
                   - c.kC1_on * S2b * E1b 
                   + c.kC1_off*C1;
            
            // E1b
            k[2] = c.kE1_on * (c.N_1 * (Pba + Pbse) / c.N_Avo 
                   - (S1b + E1b + C1 + C2))* E1 
                   - c.kE1_off * E1b
                   - c.kC1_on * S2b * E1b
                   + (c.kC1_off + c.kC1_cat)*C1
                   + c.kC2_cat*C2;
            
            // E2b
            k[3] = c.kE2_on * (c.N_2 * (Pba + Pbse) / c.N_Avo 
                   - (S2b + E2b + C1 + C2))* E2
                   - c.kE2_off * E2b
                   - c.kC2_on * S1b * E2b
                   + (c.kC2_off + c.kC2_cat)*C2
                   + c.kC1_cat * C1;

            // C1
            k[4] = c.kC1_on * S2b * E1b 
                   -(c.kC1_off + c.kC1_cat)*C1;

            // C2
            k[5] = c.kC2_on * S1b * E2b 
                   -(c.kC2_off + c.kC2_cat)*C2;
        }
};
