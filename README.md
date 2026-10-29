## MOHID Lagrangian - v26.10

MOHID Lagragian is a comprehensive high-performance Lagrangian tracer model, with sources, sinks, particle types and several options for forcing and I/O.
Altough mainly developed for oceanographic and fluvial contexts, application to atmospheric and other planetary settings should be trivial. 

Available functionalities are

- Robust pre-processing, modelling and post-processing tools
- Support for netcdf-cf files with currents, winds and wave fields, as well as water quality (salinity, temperature)
- Ability to model passive, bouyant and degrading tracers
- Stokes drift, windage, beaching, resuspension and turbulent diffusion models and options
- Ability to model millions of tracers in a modest laptop machine, paralelized processing
- Simple and fully documented simulation set-up files, ready to be abstracted by a UI
- Raw vtk time encoded output, directly compatible with Paraview and other standard post-processors and renderers
- Flexible python post processor, using cross-simulation reuseable post-processing recipes, ready to be automated
- Computation of volumetric averages and cumulative integrations, exporting the results to standard netcfd files, so you can explore the results using GIS software or publish to a thredds server
- Production of high-quality mapped plots and shapefiles using matplotlib and pandas, allowing for arbitrary calendar, integration types, subdomains including polygons and plot type combinations
- Documentation on instalation, code structure, case preparation, post processing and general useage. Fully self contained examples to get you started
- Pre-built windows executable
- Cross-platform compliant, tested and deployed
- Cmake based project, easy to set up for local compilation if required
  
## What’s new in v26.10

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

See [CHANGELOG.md](CHANGELOG.md) for the release history.

## Installation and quick start

- [Installation guide: steps 0–3](buildScripts/0-1-2-3-Installation-Guide.md)
- [Quick-start files for Linux and Windows](QuickStart/QuickStart-v.0.0/readme.md)
- [Linux installation script](MOHIDLagrangian_install.sh)

## Output examples

See the [RUN_Cases](RUN_Cases/) folder for example simulation setups and run scripts.

## Related resources

- [MOHID](http://www.mohid.com)
- [CMEMS / Copernicus Marine Service](http://marine.copernicus.eu/)
- [Code documentation](https://mohid-water-modelling-system.github.io/MOHID-Lagrangian/)

## Help, Bugs, Feedback
If you need help with MOHIDLagrangian or MOHID, want to keep up with progress, chat with developers or ask any other questions about MOHID, you can hang out by mail: <general@mohid.com> or consult our [MOHID wiki](http://wiki.mohid.com). You can also subscribe to our [MOHID forum](http://forum.mohid.com). To report bugs, please create a GitHub issue or contact any developers. More information consult <http://www.mohid.com>

## License
GNU General Public License. See the [GNU General Public License](http://www.gnu.org/copyleft/gpl.html) web page for more information.

<!--[![Build Status](https://travis-ci.org/RBCanelas/MOHID-Lagrangian.svg?branch=master)](https://travis-ci.org/RBCanelas/MOHID-Lagrangian)-->

[![License](https://img.shields.io/badge/license-GNU%20GeneraL%20Public%20License%20v3,%20GPLv3-blue.svg)]()
