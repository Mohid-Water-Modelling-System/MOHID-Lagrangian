# Changelog

## v26.10

### New features and improvements

- Added [MOHID-Lagrangian-CodeVerification](MOHID-Lagrangian-CodeVerification/README.md), with analytical and statistical verification cases, synthetic NetCDF inputs, run scripts, and Python analysis producing tables and figures.

| Verification case | Description |
|---|---|
| [CodeVerification-Advection](MOHID-Lagrangian-CodeVerification/CodeVerification-Advection/README.md) | Analytical trajectories and temporal convergence of the integration schemes. |
| [CodeVerification-Buoyancy-deriven](MOHID-Lagrangian-CodeVerification/CodeVerification-Buoyancy-deriven/README.md) | Settling and rising under constant water properties. |
| [CodeVerification-CoupledProcesses-Adv+Wd+Wv](MOHID-Lagrangian-CodeVerification/CodeVerification-CoupledProcesses-Adv+Wd+Wv/README.md) | Combined advection, windage and Stokes drift with neutral buoyancy. |
| [CodeVerification-Degradation](MOHID-Lagrangian-CodeVerification/CodeVerification-Degradation/README.md) | Analytical comparisons of detritus mass, volume and radius at constant temperature. |
| [CodeVerification-StokesDrift](MOHID-Lagrangian-CodeVerification/CodeVerification-StokesDrift/README.md) | Prescribed Stokes drift and geographical coordinate conversion. |
| [CodeVerification-Turbulence](MOHID-Lagrangian-CodeVerification/CodeVerification-Turbulence/README.md) | Random-walk and random-flight statistical checks. |
| [CodeVerification-WindDrift](MOHID-Lagrangian-CodeVerification/CodeVerification-WindDrift/README.md) | Windage scaling and geographical coordinate conversion. |
| [CodeVerification-WaterProperties](MOHID-Lagrangian-CodeVerification/CodeVerification-WaterProperties/README.md) | Independent density and dynamic-viscosity options, fallback behavior and legacy XML compatibility. |

- Refactored kernels to move shared module state into instance-owned components for improved OpenMP isolation.
- Refactored `kernelSharedProcesses` to separate common tracer processes and make them reusable across modules.
- Added separate `WindVelocity` and `StokesDriftVelocity` functions for use by `Windage` and `StokesDrift`.
- Added `DiffusionCoeffHorz` and `DiffusionCoeffVert` for independent horizontal and vertical diffusion coefficients.
- Added `UseConstantWaterDensity` with `MeanDensity` to select constant or variable seawater density.
- Added `UseConstantWaterDynamicViscosity` with `MeanDynamicViscosity` to select constant or variable seawater dynamic viscosity independently of density.
- Added `NumBlocks` to configure the number of computational domain blocks. Setting `NumBlocks` greater than or equal to the selected OpenMP thread count allows enough block tasks for all threads; performance also depends on particle distribution and workload. See [Arousa_3D_Case_Plastic_NewParam](RUN_Cases/Mohsen_Cases/Arousa_3D_Case_Plastic_NewParam/Arousa_3D_Def.xml) for configuration examples.
- Added support for reading sea surface height from hydrodynamic NetCDF files, controlled by `SshinNetcdf` through `setSshinNetcdf`.
- Added the `StokesDriftCoeff` setting for scaling the Stokes drift contribution.
- Added CPU-time reporting alongside elapsed wall time.
- Added Windows runtime dependencies and compilation scripts.
- Updated the CMake build configuration.
- Updated the [docs](docs/) documentation files and added [QuickStart](QuickStart/) files.
- Added the [installation guide for steps 0–3](buildScripts/0-1-2-3-Installation-Guide.md).
- Added [MOHIDLagrangian_install.sh](MOHIDLagrangian_install.sh).

### Bug fixes

- Corrected the corrector evaluation time in the second-order Euler predictor–corrector method.
- Fixed random-flight sampling.
- Fixed the HDF5 time-axis check for MOHID-style time arrays.
- Fixed the beaching function.
- Fixed reading NetCDF time coordinates expressed in minutes in `NetCDFparser.f90`.
- Fixed curvilinear interpolation indexing and boundary handling.
