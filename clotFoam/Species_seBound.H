class Species_seBound : public Species
{
    private:
        label injuryPatchID; 
        const scalar Kc0_on; //value only, no dimensions!
        const scalar Kc0_off; //value only, no dimensions!
        const scalar Kc0_cat; //value only, no dimensions!
        Species *fluidPhase {nullptr};

    public: 
        //- Constructor
        Species_seBound
        (
            const int& numSpecies_,
            const label& injuryPatchID_,
            const scalar Kc0_on_, 
            const scalar Kc0_off_, 
            const scalar Kc0_cat_
        ) 
        : 
            Species(numSpecies_), 
            injuryPatchID(injuryPatchID_),
            Kc0_on(Kc0_on_), Kc0_off(Kc0_off_), Kc0_cat(Kc0_cat_)
        {}

        void setPtrs(Species& fluidPhase_)
        {
            fluidPhase = &fluidPhase_;
        }

        void setInitialConcentration()
        {
            // Set the initial concentration for E0 on injuryWalls
            field[0].correctBoundaryConditions();
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
            seBound[0] -> E0,   fluidPhase[0] -> S1,
            seBound[1] -> C0,   fluidPhase[1] -> S2,
                                fluidPhase[2] -> E1,
                                fluidPhase[3] -> E2,

            input refers to the function input from RK4, 
            e.g. input[0] = E0_0 + h/2 * k1

            Output is k[i] = R_i(input)
            \* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

            const auto& E0 = input[0].boundaryField()[injuryPatchID];
            const auto& C0 = input[1].boundaryField()[injuryPatchID];
            const auto& S1 = fluidPhase->input[0].boundaryField()[injuryPatchID];

            // As far as I know, there are not dimensions for .boundaryField()
            // therefore the kinetic rates need to be passed in as scalars
            // and not dimensionedScalars.  I could be wrong, but this works.
            
            forAll(k[0].boundaryField()[injuryPatchID], faceI)
            {
                k[0].boundaryFieldRef()[injuryPatchID][faceI] =
                    -1.*Kc0_on * S1[faceI] * E0[faceI]
                    + (Kc0_off + Kc0_cat) * C0[faceI];
            
                k[1].boundaryFieldRef()[injuryPatchID][faceI] =
                    Kc0_on * S1[faceI] * E0[faceI]
                    - (Kc0_off + Kc0_cat)*C0[faceI];    
            }     
        }
};