class Species
{
    private:
        const int numSubSpecies; // number of sub-species (e.g. 4 platelet species)
    
    public: 
        //- Constructor
        Species(int numSubSpecies_) : numSubSpecies(numSubSpecies_) {}
        //- Destructor
        virtual ~Species(){}

        // public variables shared by all Species objects
        PtrList<volScalarField> field; // PtrList of sub-species solution
        PtrList<volScalarField> fieldOldTime; // PtrList of sub-species at prev time
        PtrList<volScalarField> k; //For RK solver
        PtrList<volScalarField> input; //For RK solver

        //- Member functions shared by all Species objects
        // Return the number of species for an object
        int size()
        {
            return numSubSpecies;
        }

        // Set the size of each PtrList for the object, initialize RK fields
        void setPtrListSize()
        {
            field.setSize(numSubSpecies);
            fieldOldTime.setSize(numSubSpecies);
            k.setSize(numSubSpecies);
            input.setSize(numSubSpecies);
        }

        // Initialize all fields used in the RK method to zero
        void initializeRKFields()
        {
            // scalar of dimTime for initializing k properly
            dimensionedScalar h ("h", dimTime, 1.0);

            forAll(field,j)
            {
                fieldOldTime.set(j, volScalarField(0*field[j]).clone());
                k.set(j, volScalarField(0*field[j]/h).clone());
                input.set(j, volScalarField(0*field[j]).clone());
            }
        }

        // Update the value of fieldOldTime for all species
        void setFieldOldTime()
        {
            forAll(fieldOldTime,j)
            {
                fieldOldTime[j] = field[j];
            }
        }

        // Function that edits the k's for RK4 method.  This assumes that 
        // all k's are defined as 
        //      k[i] = R_i(input) 
        // where R represents the reactions in the ADR equation 
        //      d/dt(C_i) = - div(u C_i - D grad(C_i)) + R_i(C1,C2,...)
        // Output is k[i] = R_i(input)
        virtual void updateKs() = 0;
};
