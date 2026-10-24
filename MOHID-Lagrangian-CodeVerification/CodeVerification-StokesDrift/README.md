# Stokes drift

One passive base tracer moves under a spatially uniform, constant prescribed
surface Stokes velocity (VSDX,VSDY) = (0.2,0.4) m/s. StokesDriftCoeff = 0.03
therefore gives an additive particle velocity (U,V) = (0.006,0.012) m/s.
The initial position is (1 degree, 1 degree, 0 m).

The model reads zero-current hydrodynamic NetCDF and a separate wave NetCDF
containing the prescribed drift components. Base tracers have no buoyancy
kernel. Diffusion and windage are zero, and no wind input is supplied.
Sea surface height is zero. At z = 0, the depth attenuation factor is 1.
The domain covers longitude and latitude 0.9 to 1.1 degrees, with bottom
at -200 m, leaving the full one-hour trajectory away from horizontal boundaries.

The case uses Euler with a 60 s time step and output interval. The table reports
positions at 0, 0.5 and 1 hour. The wave file also includes neutral auxiliary
wave parameters required by the input reader; the presence of both VSDX and
VSDY selects the prescribed-velocity branch. This test checks that branch,
coefficient scaling, additive transport and geographical conversion. It does
not verify the fallback calculation from wave height/period/direction or
subsurface attenuation.

## Run

Use the MOHID-Lagrangian environment and a current model executable. The
included inputs are ready. Recreate them when needed:

```bash
python analyze_results.py --prepare-inputs
```

Run and analyze separately:

```bash
bash RunCase_StokesDrift.sh
python analyze_results.py
```

The shell script recreates `StokesDrift_Def_out`. The analyzer never runs the
model. It writes `verification_table.tex`, `StokesDrift_verification.pdf`/`.png`,
`trajectory.csv` and `summary.json` inside that output folder. Use
`--input /path/to/output` for another existing output directory.

## Analytical reference

For longitude lambda and latitude phi in radians, R = 6378137 m:

```text
U = 0.03 * 0.2 = 0.006 m/s
V = 0.03 * 0.4 = 0.012 m/s
dphi/dt    = V/R
dlambda/dt = U/[R*cos(phi)]
phi(t)     = phi0 + V*t/R
lambda(t)  = lambda0 + (U/V) * [atanh(sin(phi(t))) - atanh(sin(phi0))]
z(t)       = 0
```

The reference accounts for changing latitude in the zonal metric. Angles are
converted to degrees for comparison with VTU coordinates. Long-double
arithmetic reduces cancellation in the analytical expression.
Time equals particle age: Stokes drift is active from the first Euler step.
Calendar End includes the scheduler and final-output guards.

Relative errors are normalized by analytical displacement from the initial
position. Vertical relative error is -- because analytical vertical displacement
is zero; absolute vertical motion is checked. The continuous analytical
trajectory and discrete Euler trajectory can differ slightly in longitude
because Euler evaluates the latitude-dependent metric at each step's start.
Figure displacements use the initial-latitude tangent-plane scale for display;
the analytical trajectory uses the varying latitude above.

Passing requires one persistent, wet, nonparticulate tracer, zero resolved
current and stochastic velocity, the initial and both checkpoint snapshots,
no boundary encounter, horizontal error <= 1e-3 m, vertical error <= 1e-10 m,
and checkpoint relative horizontal errors <= 1.4e-6 percent. The analyzer exits
nonzero on failure. No production Fortran changes are needed.
