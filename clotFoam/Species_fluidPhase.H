class Species_fluidPhase : public Species
{
    private:
        const chemConstants& c;
        Species *Plt {nullptr};
        Species *pltBound {nullptr};

    public: 
        //- Constructor
        Species_fluidPhase
        (
            const int& numSpecies_,
            const chemConstants& c_
        ) 
        : 
            Species(numSpecies_),
            c(c_),
            Plt(nullptr),
            pltBound(nullptr)
        {}

        void setPtrs(Species& Plt_, Species& pltBound_)
        {
            Plt = &Plt_;
            pltBound = &pltBound_;
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
            fluidPhase[0] = S1,     pltBound[0] -> S1b,     Plt[0] -> Pmu,
            fluidPhase[1] = S2,     pltBound[1] -> S2b,     Plt[1] -> Pma,
            fluidPhase[2] = E1,     pltBound[2] -> E1b,     Plt[2] -> Pba,
            fluidPhase[3] = E2,     pltBound[3] -> E2b,     Plt[3] -> Pbse.
                                    pltBound[4] -> C1,
                                    pltBound[5] -> C2,
            
            input refers to the function input from RK4, 
            e.g. input[0] = S1 + h/2 * k1

            Output is k[i] = R_i(input)
            \* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

            // fluidPhase Reference Vars
            const volScalarField& S1 = input[0];
            const volScalarField& S2 = input[1];
            const volScalarField& E1 = input[2];
            const volScalarField& E2 = input[3];

            // pltBound Reference Vars
            const volScalarField& S1b = pltBound->input[0];
            const volScalarField& S2b = pltBound->input[1];
            const volScalarField& E1b = pltBound->input[2];
            const volScalarField& E2b = pltBound->input[3];
            const volScalarField& C1 = pltBound->input[4];
            const volScalarField& C2 = pltBound->input[5];
            
            // Platelet Reference Vars
            const volScalarField& Pba = Plt->field[2];
            const volScalarField& Pbse = Plt->field[3];

            // S1
            k[0] = - c.kS1_on * (c.N_1 * (Pba + Pbse) / c.N_Avo
                   - (S1b + E1b + C1 + C2)) * S1
                   + c.kS1_off * S1b;

            // S2
            k[1] = - c.kS2_on * (c.N_2 * (Pba + Pbse) / c.N_Avo 
                   - (S2b + E2b + C1 + C2))* S2 
                   + c.kS2_off * S2b;

            // E1 (has additional reactions on SE through BC)
            k[2] = - c.kE1_on * (c.N_1 * (Pba + Pbse) / c.N_Avo 
                   - (S1b + E1b + C1 + C2))* E1 
                   + c.kE1_off * E1b;

            // E2
            k[3] = - c.kE2_on * (c.N_2 * (Pba + Pbse) / c.N_Avo 
                   - (S2b + E2b + C1 + C2))* E2
                   + c.kE2_off * E2b;
        }
};
