/**************************************************************************
    React the Chemicals M times
***************************************************************************/
{    
    //Info<< "Chemical Reaction Solver: M_react = " << M_rxn << 
    //        ", dt_react = " << h_rxn.value() << endl;

    // Solve for chemical reactions M_rxn times
    for (int m = 0; m < M_rxn; m++ )
    {
        //- RK4 Solver
        odeSolver(h_rxn, ptrToChems, num_chemObjs);
    }
}

// Update the BC's for the fluidPhase Species
forAll(fluidPhase.field, i)
{
    fluidPhase.field[i].correctBoundaryConditions();
}
