/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    clotFoam

Description
    This solver simulates blood clotting in any type of domain.
    The solver can be broken down into 3 major compents:
    1) Fluid:  Transient solver for incompressible, laminar flow of Newtonian 
       fluid with an additional Darcy term

        du/dt = - grad(p') - div[u*grad(u) - nu*grad(u)] - nu*alpha(theta_B)*u,
        div(u) = 0,

    2) Platelet Aggregation: hindered transport of 4 platelet species with 
       activation by ADP and thrombin (e2)
       
         dPmu/dt = - div[W(theta_T)*(u*Pmu - Dp*grad(Pmu))] 
                   + Rmu(Pmu,Pma,Pba,Pbse),
         dPma/dt = - div[W(theta_T)*(u*Pma - Dp*grad(Pma))] 
                   + Rma(Pmu,Pma,Pba,Pbse),
        dPba/dt  = Rba(Pmu,Pma,Pba,Pbse),
        dPbse/dt = Rbse(Pmu,Pma,Pba,Pbse),
       d[ADP]/dt = - div[u*[ADP] - Dp*grad([ADP])] + sigma_release(Pba,Pbse)

    3) Coagulation: Thrombin (E2) generation via enzymatic reactions
    
       S1 + E0 <=> C0 -> E0 + E1,
       S1 + P1 <=> S1b
       E1 + P1 <=> E1b,
       S2 + P2 <=> S2b,
       E2 + P2 <=> E2b,
       S2b + E1b <=> C1 -> E1b + E2b
       S1b + E2b <=> C2 -> E1b + E2b

    Notes: 
        - The reaction zone must be defined as a patch called "injuryWalls".
          This can be done in blockMesh or with the topoSet tool.
        - The reactive boundary conditions for the fluidPhase species are
          specified in the $FOAM_CASE/0 directory for that species. 

Author: David Montgomery 
        PhD Candidate at Colorado School of Mines 2022
        with help from Federico Municchi PhD and Karin Leiderman PhD
\*---------------------------------------------------------------------------*/

// Classes from OpenFOAM
#include "fvCFD.H"
#include "pisoControl.H"
#include "mathematicalConstants.H"

// Classes/structures for managing various species
#include "plateletConstants.H"
#include "chemConstants.H"
#include "Species_baseClass.H"
#include "Species_platelet.H"
#include "Species_seBound.H"
#include "Species_fluidPhase.H"
#include "Species_pltBound.H"

// RK4 Solver
#include "odeSolver.H"

// Rate of ADP release bell function R(tau)
double R_ADP(const double& tau)
{
    return std::exp(-1.*std::pow( tau - 3.0, 2.0)) 
        / std::sqrt(constant::mathematical::pi);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    // Create the constants and fields for simulation
    #include "createConstants.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // Set up the sigma release source term for ADP
    #include "initSigmaReleaseADP.H"

    // Set necessary pointers for each species
    #include "setSpeciesPointers.H"
    
    // Calculate initial Theta_T, Theta_B
    Plt.updateFractions();

    //--- Start time loop
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {   
        if (runTime.write())
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;
        }

        // Variable time step control variables and adjustments
        const bool adjustTimeStep =
            runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

        scalar maxCo =
            runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

        scalar maxDeltaT =
            runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);
        #include "CourantNo.H"
        #include "setDeltaT.H"

        // Solve the Navier-Stokes-Brinkman Equations
        #include "solveFluids.H"

        // Transport the platelets dp/dt = - div(W*J)
        #include "plateletTransport.H" // Transport the mobile platelets

        // Solve the reaction equations
        h_rxn = runTime.deltaT()/M_rxn; // update the reaction time-step size

        #include "plateletReactions.H"
        Plt.updateFractions();

        if (coagReactionsOn)
        {
            #include "fluidPhaseChemTransport.H"
            #include "chemReactions.H"
        }

        // Transport ADP and update sigma_release
        #include "ADP.H" 

        runTime.write();

        if (runTime.write())
        {
             // Calculate shearRate (used in shear-dependent fxns for Plt reactions)
            shearRate = Foam::sqrt(2.0) * mag(symm( fvc::grad(U) )) ;
            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;
            Info<< "max(shearRate) = "<< max(shearRate).value() <<" 1/s"<< nl << endl;
        }

        // Check if solution is diverging
        #include "isSolutionDiverging.H"
      
    } // End time loop

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
