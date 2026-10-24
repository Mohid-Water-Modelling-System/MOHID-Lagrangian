# Coupled advection, windage and Stokes drift

One paper tracer moves under simultaneous spatially uniform, constant currents,
wind-induced drift and prescribed Stokes drift. The table compares positions
at 0, 0.5 and 1 hour. Euler uses a 60 s time step and output interval.

| Contribution | Input (east, north), m/s | Coefficient | Effective velocity, m/s |
|---|---|---:|---|
| Advection | (0.20, 0.30) | 1 | (0.20, 0.30) |
| Windage | (10, 15), wind at 10 m | 0.02 | (0.20, 0.30) |
| Stokes drift | (0.2, 0.4), prescribed surface drift | 0.03 | (0.006, 0.012) |
| Total horizontal transport | | | **(0.406, 0.612)** |

The tracer starts at (1 degree, 1 degree, 0 m). The paper type activates the
buoyancy kernel, but particle and water densities are both 1026 kg/m3, so its
vertical contribution is zero. Dynamic viscosity is fixed at 0.00111834 Pa s,
giving nu = 1.09e-6 m2/s. Radius is 1 mm. Degradation, resuspension and all
diffusivities are zero. Paper does not activate plastic biofouling.

At z = 0 and zero sea surface height, wind and Stokes depth attenuation factors
are 1. Separate hydrodynamic, meteorology and wave NetCDF files are read by the
full model. The all-wet domain covers longitude and latitude 0.9 to 1.1 degrees,
with bottom at -200 m, leaving the one-hour trajectory clear of boundaries.

## Run

Use the MOHID-Lagrangian environment and a current model executable. Inputs
are included; recreate them when needed:

```bash
python analyze_results.py --prepare-inputs
```

Run and analyze separately:

```bash
bash RunCase_CoupledProcesses.sh
python analyze_results.py
```

The shell script recreates `CoupledProcesses_Def_out`. The analyzer never runs
the model. It writes `verification_table.tex`, `CoupledProcesses_verification.pdf`
and `.png`, `trajectory.csv`, and `summary.json` inside the output folder.
Use `--input /path/to/output` for another completed run. The material properties
are read from the snapshot saved with that run.

## Analytical reference

For longitude lambda and latitude phi in radians, R = 6378137 m:

```text
U = 0.20 + 0.02*10 + 0.03*0.2 = 0.406 m/s
V = 0.30 + 0.02*15 + 0.03*0.4 = 0.612 m/s
dphi/dt    = V/R
dlambda/dt = U/[R*cos(phi)]
phi(t)     = phi0 + V*t/R
lambda(t)  = lambda0 + (U/V) * [atanh(sin(phi(t))) - atanh(sin(phi0))]
z(t)       = 0
```

The analytical longitude accounts for continuously varying latitude. Positions
are converted to degrees for comparison with VTUs. Time is particle age;
horizontal transport is active from the first Euler step. Calendar End includes
scheduler and output guards. Relative errors use analytical displacement from
the initial position. Vertical relative error is -- because the analytical
vertical displacement is zero; absolute vertical motion is checked instead.

The numerical kernel adds the transport contributions before temporal integration.
This case checks their simultaneous additive composition and coordinate
conversion. It is not a test of sequential operator splitting or of every
possible process coupling. Buoyancy is evaluated at zero density contrast;
nonzero settling and rising are checked in the separate buoyancy cases.
The wave input exercises the prescribed-Stokes branch, not its wave-parameter
fallback. Diffusion and boundary interactions are outside this case's scope.

Passing requires one persistent, wet paper tracer, constant material properties,
the prescribed resolved current, no stochastic velocity or boundary encounter,
initial and both checkpoint outputs, maximum horizontal error <= 1e-3 m,
vertical error <= 1e-10 m, and checkpoint relative horizontal errors <= 1e-5
percent. The analyzer exits nonzero on failure. The expected small longitude
error comes from Euler evaluating the latitude-dependent metric at the start
of each step. No production Fortran changes are needed.
