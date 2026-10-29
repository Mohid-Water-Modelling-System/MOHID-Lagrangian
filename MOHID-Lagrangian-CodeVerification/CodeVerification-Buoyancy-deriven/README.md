# Buoyancy-driven motion

Prepared by Mohsen Shabani | Email: shabani.mohsen@outlook.com

Two cases check vertical motion under constant water properties and zero
resolved current. The existing settling case keeps its physical settings;
its name is now `Buoyancy_settling_Def`. The second case checks rising.
Both use one spherical paper tracer, which activates buoyancy without the
plastic biofouling kernel. Degradation, resuspension, diffusion, windage and
Stokes drift are disabled.

| Parameter | Settling | Rising |
|---|---:|---:|
| XML | Buoyancy_settling_Def.xml | Buoyancy_rising_Def.xml |
| Diameter | 0.002 m | 0.002 m |
| Particle density | 1050 kg/m3 | 1002 kg/m3 |
| Water density | 1026 kg/m3 | 1026 kg/m3 |
| Dynamic viscosity | 0.00111834 Pa s | 0.00111834 Pa s |
| Kinematic viscosity, mu/rho | 1.09e-6 m2/s | 1.09e-6 m2/s |
| Initial longitude, latitude | 1 degree, 1 degree | 1 degree, 1 degree |
| Initial z | 0 m | -100 m |
| Euler time step / output interval | 60 s | 60 s |
| Motion duration | 1 hour | 1 hour |

Both cases share the all-wet, zero-current NetCDF domain: longitude and latitude
0.9 to 1.1 degrees, bottom -3000 m. The rising particle stays below the surface
for the full test. Both water-property switches are 1. Constant temperature and
salinity alone would not guarantee constant density because pressure varies
with depth. The model computes the terminal velocity; it is not imposed as a
vertical current.

## Run

Use an executable built with the current water-property options. The included
NetCDF is ready. Recreate it only when needed:

```bash
python analyze_results.py --prepare-inputs
```

Run each case separately in the MOHID-Lagrangian environment:

```bash
bash RunCase_Buoyancy_settling.sh
bash RunCase_Buoyancy_rising.sh
```

Each script recreates its own `<XML name>_out` folder. After the model runs,
analyze both outputs or select one:

```bash
python analyze_results.py
python analyze_results.py --case settling
python analyze_results.py --case rising
```

For another existing output folder, use `--case rising --input /path/to/output`
(or `--case settling`). Analysis never launches the model. Each output folder
gets its own `trajectory.csv`, `summary.json`, `verification_table.tex`, and
`Buoyancy_verification.pdf`/`.png`. The two LaTeX tables have distinct labels
so they can be included in the same paper. Material properties are read from
the snapshot saved with each run.

## Analytical reference

With d the diameter, R = abs(1 - rho_p/rho_sw), and nu = mu/rho_sw:

```text
Dstar = d^3 * g * R / nu^2
speed = (nu/d) * Dstar * (38.1 + 0.93 * Dstar^(4/7))^(-7/8)
vz = sign(rho_sw - rho_p) * speed
z(t) = z0 + vz*t
```

Gravity magnitude is 9.81 m/s2. Equal absolute density contrasts give equal
speed magnitudes, with opposite directions:

| Quantity | Settling | Rising |
|---|---:|---:|
| Vertical velocity (m/s) | -0.0149925278821 | +0.0149925278821 |
| Analytical z at 0.5 h (m) | -26.986550188 | -73.013449812 |
| Analytical z at 1 h (m) | -53.973100376 | -46.026899624 |

The buoyancy kernel returns zero at particle age zero. The first Euler step
advances age but leaves position unchanged. Motion time is therefore
max(age - dt, 0), not the XML calendar clock or unadjusted age. The analyzer
checks this initial resting step. Calendar End includes it and scheduler guards
so that exact 0.5-hour and 1-hour motion snapshots are available. Keep dt,
OutputWriteTime, source scheduling and End consistent if edited.

Vertical relative error is normalized by analytical vertical displacement.
Horizontal relative errors are -- because the analytical displacement is zero;
absolute horizontal drift is checked instead.

Passing requires one persistent particle, zero resolved velocity, a wet state,
constant material properties, no boundary encounter, both checkpoint outputs,
maximum vertical error <= 1e-6 m, and horizontal drift <= 1e-10 degrees.
The script exits nonzero on failure. These cases verify constant-property
terminal motion and vertical integration. Variable-property motion is separate.
