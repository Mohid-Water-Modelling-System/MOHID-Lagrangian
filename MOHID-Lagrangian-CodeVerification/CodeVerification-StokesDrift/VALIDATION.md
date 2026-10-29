# Stokes-drift verification

Prepared by Mohsen Shabani | Email: shabani.mohsen@outlook.com

The preprocessor and full one-hour run completed. The VTU analyzer returned
PASS. No production Fortran equations were changed.

| Time (h) | Analytical longitude | Numerical longitude | Analytical latitude | Numerical latitude | Longitude error (%) | Latitude error (%) |
|---|---:|---:|---:|---:|---:|---:|
| 0.5 | 1.0000970328321 | 1.0000970328320 | 1.0001940361014 | 1.0001940361014 | 9.8399e-08 | 1.6021e-09 |
| 1 | 1.0001940656700 | 1.0001940656698 | 1.0003880722027 | 1.0003880722027 | 9.8513e-08 | 1.6593e-09 |

Maximum horizontal error: 2.12909e-08 m.
Vertical position remains exactly 0 m. No boundary encounter, nonzero resolved
current or stochastic velocity was detected.

All hydrodynamic u/v/w values are exactly zero. All wave VSDX/VSDY values
are exactly 0.2/0.4 m/s. StokesDriftCoeff = 0.03 gives the effective particle
velocity 0.006/0.012 m/s. Surface depth attenuation is 1.

The analytical longitude was independently checked against 32-point Gaussian
quadrature within 2e-15 degrees at both checkpoints. The small longitude
difference is consistent with Euler evaluating the latitude-dependent metric
at the start of each 60 s step. Latitude has constant angular velocity; its
remaining discrepancy is consistent with floating-point roundoff.

This verifies the prescribed-component Stokes branch and its scaled additive
contribution. It does not test the wave-parameter fallback or depth attenuation.

The LaTeX table, PDF/PNG figure, CSV and JSON summary are in
`StokesDrift_Def_out`. The shell script replaces that folder on the next run.
