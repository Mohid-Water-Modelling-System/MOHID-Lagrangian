# Buoyancy verification

Prepared by Mohsen Shabani | Email: shabani.mohsen@outlook.com

Both one-hour cases pass the trajectory checks. Settling uses the preserved
output from the earlier run, renamed to `Buoyancy_settling_Def_out`.
Rising was run through the preprocessor and full model with the new XML.
No production Fortran equations or settling-case physical settings were changed.

| Case | Motion time (h) | Analytical z (m) | Numerical z (m) | Relative error (%) |
|---|---:|---:|---:|---:|
| settling | 0.5 | -26.986550187777 | -26.986550187777 | 2.633e-14 |
| settling | 1 | -53.973100375554 | -53.973100375554 | 2.633e-14 |
| rising | 0.5 | -73.013449812223 | -73.013449812223 | 5.2659e-14 |
| rising | 1 | -46.026899624446 | -46.026899624446 | 3.9494e-14 |

Maximum vertical errors: 1.4211e-14 m for settling and 2.1316e-14 m for rising.
Horizontal drift is zero in both cases. Both particles remain wet, away from
the relevant boundary, with constant material properties.

Motion time excludes the first age-zero resting step. Checkpoints at 0.5 and
1 hour have particle ages 1860 and 3660 seconds for dt = 60 s.
The saved settling XML still records the original longer run; only its first
hour is analyzed. Its VTU data were preserved.

Each output folder contains its own LaTeX table, figure, CSV and JSON summary.
The water-property regression suite also passes: 80 property checks,
20 invalid-setting checks and 7 XML compatibility checks.
