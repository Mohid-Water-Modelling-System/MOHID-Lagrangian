# Advection verification

Prepared by Mohsen Shabani | Email: shabani.mohsen@outlook.com

One passive `base` particle starts at `(X,Y,Z) = (100,100,-100)` m in a
2 km × 2 km water domain, 400 m deep. Only current-driven advection is active.
All diffusion coefficients and wind/wave contributions are zero.

The shared NetCDF contains

```
u = a X cos(latitude) (1 + gamma t)
v = b Y (1 + gamma t)
w = -(a+b) Z (1 + gamma t)
```

with `a=5e-5 /s`, `b=2.5e-5 /s`, `gamma=1e-3 /s` and physical duration
3600 s. The longitude metric factor cancels the geographical conversion in
the model, so the local coordinate trajectories are exponential. Here
`X=degree*longitude`, `Y=degree*latitude`, `Z=depth`, with the same Earth
radius and initialized pi as `utilities.f90`. This is a local coordinate
verification, not an assertion of Cartesian incompressibility on the sphere.

The grid has 21 longitudes, 201 latitudes and nine levels, with forcing every
60 s. Linear interpolation exactly reproduces the linear factors; the
latitude grid resolves the smooth metric factor. The filename is shared by
all 12 cases. No production Fortran code is changed.

## Run

With the model compiled and Conda available, for example from this folder:

```bash
./RunCase_Advetion_Verification_Eluer1_030.sh
```

Scripts locate the repository from their own path, so they can also be
launched from another directory. Outputs remain inside the case folder.
The other scripts have the same structure. Methods are `Eluer1`, `Eluer2` and
`RK4`; time steps are `015`, `030`, `060` and `120` seconds. Names retain the
requested spelling. Each script creates its own `${name}_out` folder and
clears that folder when run again. The `data` folder follows the turbulence
case structure; only age and particulate are needed as optional output fields.

The filename defines the method and time step. To synchronize Start, End,
Integrator, output interval and source timing in all 12 source XMLs:

```bash
python analyze_results.py --sync-cases
```

This preserves XML comments and does not change existing `_out` folders.
Each RunCase script checks its XML before clearing previous results. If a
setting is inconsistent, it stops with a repair command. After synchronizing,
rerun the affected cases before analyzing their results. Do not repair the
XML inside an existing output folder: it describes the run that produced it.

If needed, regenerate the shared input before running:

```bash
python analyze_results.py --prepare-inputs
```

The initial loop does not advance the solver. Start dates are offset so the
first actual solver evaluation is at 00:04:00 on 1 January 2024 in every case.
The End date includes the solver's strict final-step guard. Each particle
therefore undergoes 3600 s of motion with the same physical forcing clock.
The initial emission stack and half-emission increment generate exactly one
particle during the only active source step.

### Timing hint

Keep the existing Start and End settings for each time step. The model skips
integration in its initial cycle and uses a strict final-step condition;
these calendar dates give every case the same 3600 s of particle motion.
All dates below are on 1 January 2024.

| Time step (s) | Start | End |
|---:|---|---|
| 15 | 00:03:45 | 01:04:15 |
| 30 | 00:03:30 | 01:04:30 |
| 60 | 00:03:00 | 01:05:00 |
| 120 | 00:02:00 | 01:06:00 |

`OutputWriteTime` controls saving frequency, not the forcing clock. Changing
it does not correct inconsistent Start/End settings. The present analysis
uses output at every integration step; keep the supplied settings.


## Analyze

After running all 12 cases:

```bash
python analyze_results.py
```

The script checks particle identity, age, water state, boundary clearance,
zero diffusion settings and the forcing clock. Each output folder receives
`trajectory.csv`, `summary.json` and `trajectory_table.tex`. The trajectory
LaTeX table compares analytical and numerical longitude, latitude and depth
at particle ages 0, 1800 and 3600 s, with component-wise percentage errors
normalized by analytical displacement. It requires `booktabs` and `graphicx`.
This folder receives:

- `convergence.csv`: signed coordinate errors, relative coordinate errors,
  Euclidean position error and observed order.
- `convergence_table.tex`: the manuscript table (requires `booktabs`, `multirow` and
  `amsmath`). The caption references the error-norm and observed-order
  equation labels in the paper.
- `Advection_convergence.pdf` and `.png`: final position error versus time step.

Relative coordinate errors use analytical displacement as the denominator.
Observed order is `log2(E_dt/E_dt/2)` on the coarser row. No expected order is
forced onto measured results. One particle suffices for this deterministic
trajectory test; it does not validate every possible initial position.

For an interpolation check, retain a copy of the original outputs, prepare
input with `--prepare-inputs --latitude-refinement 2`, and rerun selected
cases. Restore the default input afterwards. The same forcing equations and
time sampling apply to both grids.

Measured errors and the grid check are in [VALIDATION.md](VALIDATION.md).
