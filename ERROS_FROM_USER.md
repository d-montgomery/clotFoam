# Error 1 
A lot of "parse error while scanning" 

```

================================================

OpenFOAM Environment Check

================================================

WM_PROJECT_VERSION: v2412

FOAM_USER_APPBIN: /home/ofuser/OpenFOAM/ofuser-v2412/platforms/linuxARM64GccDPInt32Opt/bin

================================================


================================================

Compiling clotFoam solver...

================================================

Making dependencies: clotFoam.C

wmkdepend: parse error while scanning 'plateletConstants.H' ... perhaps missing a final newline

wmkdepend: parse error while scanning 'chemConstants.H' ... perhaps missing a final newline

wmkdepend: parse error while scanning 'Species_seBound.H' ... perhaps missing a final newline

wmkdepend: parse error while scanning 'Species_fluidPhase.H' ... perhaps missing a final newline

wmkdepend: parse error while scanning 'Species_pltBound.H' ... perhaps missing a final newline

wmkdepend: parse error while scanning 'odeSolver.H' ... perhaps missing a final newline

wmkdepend: parse error while scanning 'createFields.H' ... perhaps missing a final newline

wmkdepend: parse error while scanning 'setSpeciesPointers.H' ... perhaps missing a final newline

wmkdepend: parse error while scanning 'solveFluids.H' ... perhaps missing a final newline

wmkdepend: parse error while scanning 'chemReactions.H' ... perhaps missing a final newline

wmkdepend: parse error while scanning 'isSolutionDiverging.H' ... perhaps missing a final newline

g++ -std=c++17 -pthread -DOPENFOAM=2412 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -Wno-attributes -Wno-unknown-pragmas -O3 -floop-optimize -falign-loops -falign-labels -falign-functions -falign-jumps  -DNoRepository -ftemplate-depth-100  -I/usr/lib/openfoam/openfoam2412/src/finiteVolume/lnInclude -I/usr/lib/openfoam/openfoam2412/src/meshTools/lnInclude -iquote. -IlnInclude -I/usr/lib/openfoam/openfoam2412/src/OpenFOAM/lnInclude -I/usr/lib/openfoam/openfoam2412/src/OSspecific/POSIX/lnInclude   -fPIC -c clotFoam.C -o Make/linuxARM64GccDPInt32Opt/clotFoam.o

```

# Warning 1

```
/usr/lib/openfoam/openfoam2412/src/OpenFOAM/lnInclude/dimensionedType.C:332:1: note: declared here

  332 | Foam::dimensioned<Type>::dimensioned

      | ^~~~

plateletConstants.H:30:5: warning: 'Foam::dimensioned<Type>::dimensioned(const Foam::word&, const Foam::dimensionSet&, Foam::Istream&) [with Type = double]' is deprecated: Since 2018-11; use "construct from dictionary or entry" [-Wdeprecated-declarations]

   30 |     };

      |     ^

/usr/lib/openfoam/openfoam2412/src/OpenFOAM/lnInclude/dimensionedType.C:332:1: note: declared here

  332 | Foam::dimensioned<Type>::dimensioned

      | ^~~~

```