#!/bin/bash
set -e

echo "================================================"
echo "OpenFOAM Environment Check"
echo "================================================"
echo "WM_PROJECT_VERSION: $WM_PROJECT_VERSION"
echo "FOAM_USER_APPBIN: $FOAM_USER_APPBIN"
echo "================================================"

# Compile clotFoam
echo ""
echo "================================================"
echo "Compiling clotFoam solver..."
echo "================================================"
cd /home/ofuser/clotFoam
wclean
wmake

# Verify the executable was created
if [ ! -f "$FOAM_USER_APPBIN/clotFoam" ]; then
    echo "ERROR: clotFoam executable not found at $FOAM_USER_APPBIN/clotFoam"
    exit 1
fi

echo ""
echo "================================================"
echo "Compilation successful!"
echo "================================================"

# Run rectangle2D tutorial
echo ""
echo "================================================"
echo "Running rectangle2D tutorial..."
echo "================================================"
cd /home/ofuser/tutorials/rectangle2D

# Clean any previous results
rm -rf [1-9]* 0.* processor* log* postProcessing 2>/dev/null || true

# Run blockMesh
echo "Running blockMesh..."
blockMesh > log.blockMesh 2>&1
if [ $? -ne 0 ]; then
    echo "ERROR: blockMesh failed"
    cat log.blockMesh
    exit 1
fi

# Modify controlDict to run 5-10 timesteps for testing
cp system/controlDict system/controlDict.orig
cat > system/controlDict << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}

application     clotFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         0.0015;
deltaT          0.00005;
writeControl    timeStep;
writeInterval   5;
purgeWrite      0;
writeFormat     ascii;
writePrecision  6;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable false;
adjustTimeStep  no;
maxCo           0.25;
maxDeltaT       0.0001;
EOF

# Run clotFoam - should complete 10 timesteps quickly
# Timeout set to 5 minutes as safety margin
echo "Running clotFoam (10 timesteps for testing)..."
timeout 300 clotFoam > log.clotFoam 2>&1
CLOTFOAM_EXIT=$?

# Restore original controlDict
mv system/controlDict.orig system/controlDict

echo ""
echo "clotFoam finished with exit code: $CLOTFOAM_EXIT"

# Check if clotFoam at least started successfully
if ! grep -q "Create mesh for time" log.clotFoam; then
    echo "ERROR: clotFoam failed to initialize"
    echo "First 100 lines of log:"
    head -100 log.clotFoam
    exit 1
fi

# Check exit code
if [ $CLOTFOAM_EXIT -eq 124 ]; then
    echo "ERROR: Solver timed out - should complete 10 steps quickly"
    echo "Last 50 lines of log:"
    tail -50 log.clotFoam
    exit 1
elif [ $CLOTFOAM_EXIT -ne 0 ]; then
    echo "ERROR: clotFoam failed with exit code $CLOTFOAM_EXIT"
    echo "Last 50 lines of log:"
    tail -50 log.clotFoam
    exit 1
fi

# Check if solver completed successfully
if ! grep -q "^End" log.clotFoam; then
    echo "ERROR: Solver did not complete successfully (no 'End' marker)"
    tail -50 log.clotFoam
    exit 1
fi

# Count timesteps completed
TIMESTEPS=$(grep -c "^Time = " log.clotFoam || echo "0")
echo "✓ Completed $TIMESTEPS time steps"

if [ $TIMESTEPS -lt 5 ]; then
    echo "ERROR: Expected at least 5 timesteps, got $TIMESTEPS"
    exit 1
fi

# Verify output directory was created
if [ -d "0.001" ]; then
    echo "✓ Found output directory: 0.001"
elif [ -d "0.001" ] || [ -d "1e-03" ] || ls -d [0-9]* 2>/dev/null | grep -v "^0$" | head -1 > /dev/null; then
    LATEST_TIME=$(ls -d [0-9]* 2>/dev/null | grep -v "^0$" | sort -n | tail -1)
    echo "✓ Found output directory: $LATEST_TIME"
else
    echo "WARNING: No output time directory found, but solver completed"
fi

echo ""
echo "================================================"
echo "All tests passed successfully!"
echo "================================================"
