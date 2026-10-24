# Wind-induced drift

One passive base tracer moves under a spatially uniform, constant wind at
10 m height: u10 = 10 m/s, v10 = 15 m/s. WindDragCoeff = 0.02 gives surface
windage velocities U = 0.2 m/s and V = 0.3 m/s. The tracer starts at
(1 degree, 1 degree, 0 m).

Zero-current hydrodynamic NetCDF and constant-wind meteorology NetCDF are
read by the full model. The base tracer has no buoyancy kernel. Diffusion and
Stokes drift are zero, and no wave input is supplied. At z = 0 with zero sea
surface height, the wind depth-attenuation factor is 1.

The all-wet domain covers longitude and latitude 0.9 to 1.1 degrees and depths
-200 to 0 m. This leaves the complete one-hour trajectory far from the
horizontal boundaries. Euler and OutputWriteTime both use 60 seconds.
The table compares positions at 0, 0.5 and 1 hour.

## Run

Use the MOHID-Lagrangian environment and an executable built from the current
source. The included NetCDF files are ready. Recreate them when needed:

```bash
python analyze_results.py --prepare-inputs
```

Run the model, then analyze the output separately:

```bash
bash RunCase_WindDrift.sh
python analyze_results.py
```

The shell script recreates `WindDrift_Def_out`. Analysis never runs the model.
It writes these files inside that output folder:

- `verification_table.tex`: analytical/numerical positions and relative errors.
- `WindDrift_verification.pdf` and `.png`: displacement and position error.
- `trajectory.csv`: every output time.
- `summary.json`: physical parameters, acceptance checks and checkpoints.

Use `--input /path/to/output` to analyze another existing output folder.

## Analytical reference

Let longitude lambda and latitude phi be in radians, R = 6378137 m, U = 0.2 m/s,
and V = 0.3 m/s. The model's spherical-coordinate equations are:

```text
dphi/dt    = V/R
dlambda/dt = U/[R*cos(phi)]
phi(t)     = phi0 + V*t/R
lambda(t)  = lambda0 + (U/V) * [atanh(sin(phi(t))) - atanh(sin(phi0))]
z(t)       = 0
```

Positions are converted to degrees for comparison with VTU coordinates.
The changing latitude is included in the analytical longitude. Holding latitude
fixed would omit the geographical metric effect being checked. Long-double
arithmetic reduces cancellation when evaluating the reference expression.

Time is particle age. Windage is active from the first Euler integration step.
Calendar End includes the initial scheduler cycle and final output guard.
The analytical trajectory is continuous in time; Euler's small temporal error
in longitude is expected and is not removed from the reported comparison.
This test checks the windage implementation and coordinate conversion, not
higher-order temporal convergence.

Relative errors are normalized by analytical displacement from the release
point. Vertical relative error is -- because the analytical displacement is
zero; absolute vertical motion is checked instead. Figure displacements use
a fixed tangent-plane scale at the release latitude for display only; the
reference trajectory uses the varying latitude above.

Passing requires one persistent, wet, nonparticulate particle, zero resolved
current and stochastic velocity, initial and exact 0.5/1-hour outputs, no boundary
encounter, maximum horizontal error <= 1e-3 m, vertical error <= 1e-10 m, and
checkpoint relative horizontal errors <= 6e-6 percent. Analysis exits nonzero
on failure. No production Fortran changes are needed.
