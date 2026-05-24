# Migration from OpenFOAM 9 to OpenFOAM 12

## Overview
This document tracks all changes made to migrate clotFoam from OpenFOAM 9 to OpenFOAM 12.

## Migration Date
Started: January 14, 2026

## Breaking Changes

### API Changes

1. **PtrList::set() Overload Ambiguity**
   - **Issue**: Multiple overloads for `set()` cause ambiguity with temporary objects
   - **Solution**: Use `.clone()` on temporary field objects
   - **Example**: `fieldOldTime.set(j, volScalarField(0*field[j]).clone())`

2. **Dictionary Lookup Methods**
   - **Issue**: `lookup<T>()` method removed
   - **Solution**: Use `get<T>()` instead
   - **Example**: `runTime.controlDict().get<scalar>("startTime")`

3. **writeEntry() Function Signature**
   - **Issue**: Global `writeEntry(os, key, value)` deprecated
   - **Solution**: Use member function `value.writeEntry(key, os)`
   - **Example**: `fld.writeEntry("value", os)`

### Build System Changes

No breaking changes in build system - existing `Make/files` and `Make/options` work without modification.

### Runtime Changes

1. **Security Restrictions**
   - **Issue**: OpenFOAM 12 prevents `#calc` and `#codeStream` execution as root user
   - **Solution**: Run solver as non-root user
   - **Impact**: Docker containers must create and use non-root user

2. **Constant Naming**
   - **Issue**: `small` constant no longer available in global namespace
   - **Solution**: Use `SMALL` (uppercase) constant
   - **Example**: `maxCo/(CoNum + SMALL)`

## Change Log

### Phase 1: Docker Testing Infrastructure
- ✅ Created `docker/Dockerfile` using `opencfd/openfoam-dev:2412`
- ✅ Created `docker/test.sh` for automated testing
- ✅ Created `Makefile` with Docker convenience targets
- ✅ Created `docker-compose.yml` for local testing
- ✅ Created `.github/workflows/openfoam12.yml` for CI/CD

### Phase 2: Compilation Fixes
- ✅ Fixed ambiguous `PtrList::set()` calls in `Species_baseClass.H` - Added `.clone()` to resolve overload ambiguity
- ✅ Changed `lookup<scalar>()` to `get<scalar>()` in `initSigmaReleaseADP.H` - API change in dictionary access
- ✅ Changed `small` constant to `SMALL` in `setDeltaT.H` - Constant naming convention change
- ✅ Fixed deprecated `dimensioned` constructors in `plateletConstants.H` - Changed from `lookup()` to `get<scalar>()`
- ✅ Fixed deprecated `dimensioned` constructors in `chemConstants.H` - Changed from `lookup()` to `get<scalar>()`
- ✅ Fixed deprecated `dimensioned` constructors in `createConstants.H` - Changed from `lookup()` to `get<scalar>()`
- ✅ Added missing final newlines to 11 header files - Fixed wmkdepend parse warnings

### Phase 3: Runtime Fixes
- ✅ Fixed OpenFOAM security restriction - Created non-root user (`ofuser`) in Docker to allow `#calc` and `#codeStream` execution
- ✅ Fixed `writeEntry()` API change in tutorial files - Changed from `writeEntry(os, "", fld)` to `fld.writeEntry("value", os)` in:
  - `tutorials/rectangle2D/0/Hadh`
  - `tutorials/Tjunction2D/0/Hadh`
  - `tutorials/Hjunction3D/0/Hadh`

## Testing Strategy
- Build solver in OpenFOAM 12 Docker container
- Run rectangle2D tutorial with limited timesteps (6-10 steps, endTime=0.0015)
- Verify compilation and runtime stability
- Focus on stability, not result validation
- Tests complete in ~2-3 minutes

## Known Issues
None - All compilation and runtime issues resolved ✅

## References
- [OpenFOAM Wiki: Porting from older versions to 1.3](https://openfoamwiki.net/index.php/HowTo_port_from_older_versions_to_13)
- OpenFOAM 12 Release Notes
- clotFoam Original Paper: https://doi.org/10.1016/j.softx.2023.101483
