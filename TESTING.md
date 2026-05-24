# clotFoam Testing Guide

## Docker-Based Testing

The project includes Docker-based testing infrastructure for OpenFOAM 12 that ensures reproducible builds and test execution.

### Prerequisites

- Docker installed and running
- Docker Desktop or Docker Engine
- At least 4GB RAM allocated to Docker

### Quick Start

Test compilation and runtime:

```bash
make docker-test
```

This will:
1. Build the Docker image with OpenFOAM 12
2. Compile the clotFoam solver
3. Run the rectangle2D tutorial with limited timesteps
4. Verify successful execution

### Manual Testing

#### Build Docker Image

```bash
make docker-build
```

#### Interactive Shell

Open an interactive shell in the Docker container:

```bash
make docker-shell
```

Inside the container:

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
cd clotFoam
wmake
cd ../tutorials/rectangle2D
blockMesh
clotFoam
```

### CI/CD Integration

The project includes GitHub Actions workflow (`.github/workflows/openfoam12.yml`) that:
- Triggers on push to main/master/develop branches
- Triggers on pull requests
- Builds the solver in OpenFOAM 12
- Runs automated tests
- Archives test logs as artifacts

### Test Configuration

Test parameters are configured in `docker/test.sh`:
- **End Time**: 0.05s (limited for fast testing)
- **Time Step**: 0.0001s
- **Max Courant Number**: 0.25
- **Test Case**: rectangle2D

### Troubleshooting

#### Docker Build Fails

```bash
# Clean Docker cache
docker system prune -a

# Rebuild from scratch
docker build --no-cache -t clotfoam:openfoam12 -f docker/Dockerfile .
```

#### Compilation Errors

Check that all code changes from OF9→OF12 migration are applied:
- `Species_baseClass.H`: `.clone()` on temporary fields
- `initSigmaReleaseADP.H`: `get<scalar>()` instead of `lookup<scalar>()`
- `setDeltaT.H`: `SMALL` instead of `small`
- Tutorial `Hadh` files: `fld.writeEntry("value", os)` syntax

#### Runtime Errors

If clotFoam fails during execution:

```bash
# Check logs inside container
docker run --rm -it clotfoam:openfoam12 /bin/bash
cd /home/ofuser/tutorials/rectangle2D
cat log.clotFoam
```

### Test Logs

Test execution creates the following logs:
- `log.blockMesh`: Mesh generation output
- `log.clotFoam`: Solver execution output

These logs are available as artifacts in GitHub Actions runs.

### Performance Notes

- Docker testing uses ARM64 architecture on Apple Silicon
- Compilation takes ~30-60 seconds
- rectangle2D test case runs in ~10-20 seconds
- Full tutorial simulations can take hours (not included in automated tests)

### Advanced Usage

#### Custom Test Duration

Modify `docker/test.sh` to change test duration:

```bash
endTime         0.1;    # Run longer
maxDeltaT       0.002;  # Larger timesteps
```

#### Test All Tutorials

```bash
docker run --rm -it clotfoam:openfoam12 /bin/bash
# Manually run Tjunction2D or Hjunction3D
```

#### Parallel Testing

```bash
# Use docker-compose for parallel execution
docker-compose up
```

## Local Testing (without Docker)

Requires OpenFOAM 12 installed locally:

```bash
cd clotFoam
wclean
wmake

cd ../tutorials/rectangle2D
blockMesh
clotFoam
```

See [README.md](README.md) for detailed installation instructions.
