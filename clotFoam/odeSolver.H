/*****************************************************************************
    RK4 Method for a species object y
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
y = y0 + h/6 * (k1 + 2*k2 + 2*k3 + k4)
k0 = 0*y0,
k1 = y0 + 0.0*h*k0,
k2 = y0 + 0.5*h*k1,
k3 = y0 + 0.5*h*k2,
k4 = y0 + 1.0*h*k3 
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
This solver iterates 4 times, updating each species object y[i] in each 
iteration using newly computed k's that are scalled by the "a" coefficeint.  
This avoids storing k1, k2, k3, and k4 for each subspecies of y[i]. That is, 
the method creates one k of size y[i]->size(), instead of 4 k's each of 
size y[i]->size().  
******************************************************************************/

void odeSolver
(
    const dimensionedScalar& h, // time step size
    Species* y[], // Array of ptrs to objects of type Species
    int numSpecies // number of species objects
)
{

    // Coefficients for RK4 method. 
    double alpha[4] = {1.0, 2.0, 2.0, 1.0};
    double beta[4] = {0.0, 0.5, 0.5, 1.0};
    
    // Set fieldOldTime = currentField for each species object
    for(int i = 0; i < numSpecies; i++)
    {
        y[i]->setFieldOldTime();
    }

    // Iterate to calculate k1,...,k4 and update y
    for (int n = 0; n < 4; n++) 
    {
        // For each species object y[i] calculate: 
        // input = fieldOldTime + h * beta[n] *  k[n-1]
        // for the jth subspecies of that object
        for(int i = 0; i < numSpecies; i++)
        {
            forAll(y[i]->input, j)
            {
                y[i]->input[j] = y[i]->fieldOldTime[j] 
                                + h * beta[n] * y[i]->k[j];
            }
        }

        // For each species object get k_n = RHS(species.input)
        for(int i = 0; i < numSpecies; i++)
        {
            y[i]->updateKs();
        }

        // For each species object, i
        for(int i = 0; i < numSpecies; i++)
        {
            // Update each field[j] for ith Species object using RK4 Method
            forAll(y[i]->field, j)
            {
                y[i]->field[j] = y[i]->field[j] 
                                + h / 6. * alpha[n] * y[i]->k[j]; 
            }
        }
    }
    
}
