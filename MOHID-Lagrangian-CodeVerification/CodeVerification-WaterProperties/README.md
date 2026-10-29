# Water properties

Prepared by Mohsen Shabani | Email: shabani.mohsen@outlook.com

Select water density and dynamic viscosity independently in the XML constants:

```xml
<UseConstantWaterDensity value="1" />
<MeanDensity value="1027" units_comment="kg/m3" />
<UseConstantWaterDynamicViscosity value="0" />
<MeanDynamicViscosity value="0.001" units_comment="Pa s" />
```

Both switches default to 0. This example uses constant density and calculates
mu from temperature and salinity. MeanDynamicViscosity defaults to 0.001 Pa s.
It replaces the previous WaterDynamicViscosity XML name.

| Density switch | Dynamic viscosity switch | Density | Dynamic viscosity |
|---:|---:|---|---|
| 0 | 0 | rho(T,S,z) | mu(T,S) |
| 1 | 0 | MeanDensity | mu(T,S) |
| 0 | 1 | rho(T,S,z) | MeanDynamicViscosity |
| 1 | 1 | MeanDensity | MeanDynamicViscosity |

Each water property is checked inside its own shared function before division. Nonfinite density,
or abs(1 - rho/MeanDensity) >= 0.9, is replaced with MeanDensity. Nonfinite dynamic
viscosity, or abs(1 - mu/MeanDynamicViscosity) >= 0.9, is replaced with
MeanDynamicViscosity. The reference inputs must be finite and positive.

Kinematic viscosity is always nu = mu/rho using those corrected properties.
There is no separate kinematic-viscosity limiter. MeanKinematicViscosity remains readable
for compatibility with older XML files but no longer controls these water-property
calculations. Its units are m2/s, so it must not be copied into MeanDynamicViscosity.
With constant mu and variable density, nu varies with depth.

Missing temperature uses 15 degrees C; missing salinity uses 36 (PSS-78).
Each available field is retained independently. Calculated density includes
pressure dependence on particle depth, even when T and S are fixed defaults.
Particle density and the existing particle/water contrast safeguards are unchanged.

## Check

Compile the main project, then run:

```bash
python run_tests.py
```

The runner links the actual model library using the configured CMake compiler.
Use `--build /path/to/build` for another build directory. No particle simulations
run. Checks cover all four switches, the original equations, depth dependence,
missing fields, dynamic-viscosity guard thresholds, nonfinite inputs, density
fallback, buoyancy, XML reading and invalid settings.
`legacy_water_reference.f90` retains the original equations for comparison.
Failed checks return a nonzero exit status.

## Kinematic-viscosity input names

The internal Fortran name is `MeanKinematicViscosity`, with setter
`setMeanKinematicViscosity`. XML accepts `MeanKinematicViscosity` (preferred),
`MeanKViscosity`, and `MeanKVisco`. The preferred name takes priority when
several are present; otherwise `MeanKViscosity` takes priority over `MeanKVisco`.
The value remains a legacy input; water kinematic viscosity is calculated as mu/rho.
