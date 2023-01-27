/**************************************************************************
    Set necessary pointers for each species object
***************************************************************************/
// Set necessary pointers for each species
Plt.setPtrs(fluidPhase); //include Ptr to fluidPhase in Plt obj
seBound.setPtrs(fluidPhase); //include Ptr to fluidPhase in seBound obj
pltBound.setPtrs(Plt, fluidPhase); //Ptr to Plt and fluidPhase objs
fluidPhase.setPtrs(Plt, pltBound); //Ptr to Plt and pltBound objs

// Array of pointers to Plt object (for RK method)
int num_pltObjs = 1;
Species *ptrToPlt[num_pltObjs];
ptrToPlt[0] = &Plt;

// Array of pointers to all chemical species objects(for RK method)
int num_chemObjs = 3;
Species *ptrToChems[num_chemObjs]; // Make a list of ptrs for chem species RK solver 
ptrToChems[0] = &seBound;
ptrToChems[1] = &fluidPhase;
ptrToChems[2] = &pltBound;