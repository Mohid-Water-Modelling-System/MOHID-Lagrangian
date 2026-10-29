# Wind-drift verification

Prepared by Mohsen Shabani | Email: shabani.mohsen@outlook.com

The preprocessor and full one-hour model run completed. The VTU analyzer
returned PASS. No production Fortran changes were made.

| Time (h) | Analytical longitude | Numerical longitude | Analytical latitude | Numerical latitude | Longitude error (%) | Latitude error (%) |
|---|---:|---:|---:|---:|---:|---:|
| 0.5 | 1.0032344300356 | 1.0032344299558 | 1.0048509025342 | 1.0048509025342 | 2.469e-06 | 9.1548e-12 |
| 1 | 1.0064688648744 | 1.0064688647143 | 1.0097018050685 | 1.0097018050685 | 2.4749e-06 | 9.1548e-12 |

Maximum horizontal error: 1.78194e-05 m.
Vertical position remains exactly 0 m. The particle remained wet, with zero
resolved current and stochastic velocity, throughout the run.

The small longitude error is consistent with forward Euler evaluating the
latitude-dependent metric at the start of each step. The analytical reference
includes continuous variation of latitude. It was independently checked against
32-point Gaussian quadrature to within 2e-15 degrees at both checkpoints.

All hydrodynamic NetCDF u/v/w values were checked to be exactly zero.
All meteorology u10/v10 values were checked to be exactly 10/15 m/s.

The table, PDF/PNG figure, trajectory CSV and JSON summary are in
`WindDrift_Def_out`. The run script replaces this folder on the next execution.

## Interpretation of coordinate errors

Paper-ready LaTeX:

```latex
The difference between the longitudinal and latitudinal errors follows
from the spherical-coordinate formulation and the temporal integration
scheme. Under the prescribed constant wind, the meridional angular velocity,
$\mathrm{d}\phi/\mathrm{d}t=V/R$, is constant; forward Euler therefore
introduces no temporal truncation error in latitude, and the residual
relative error of approximately $9.15\times10^{-12}\%$ is consistent with
floating-point roundoff. In contrast, the zonal angular velocity,
$\mathrm{d}\lambda/\mathrm{d}t=U/[R\cos\phi(t)]$, varies as the particle
moves northward. Evaluating this latitude-dependent metric at the beginning
of each $60~\mathrm{s}$ step produces a longitudinal relative error of
approximately $2.47\times10^{-6}\%$ after $1~\mathrm{h}$, corresponding to
an absolute horizontal position error of $1.78\times10^{-5}~\mathrm{m}$.
Here, $U=C^{\mathrm{wind}}u_{10}$ and $V=C^{\mathrm{wind}}v_{10}$, and
$R$ is the Earth radius used by the model. The error disparity is thus
consistent with the expected temporal discretization of the spherical
metric under constant physical wind forcing.
```
