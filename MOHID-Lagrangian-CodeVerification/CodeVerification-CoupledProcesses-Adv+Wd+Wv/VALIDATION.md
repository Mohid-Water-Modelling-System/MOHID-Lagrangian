# Coupled-process verification

The preprocessor and complete one-hour model run finished. VTU analysis returned
PASS. No production Fortran equations or other verification cases were changed.

| Time (h) | Analytical longitude | Numerical longitude | Analytical latitude | Numerical latitude | Longitude error (%) | Latitude error (%) |
|---|---:|---:|---:|---:|---:|---:|
| 0.5 | 1.0065658980428 | 1.0065658977113 | 1.0098958411699 | 1.0098958411699 | 5.0492e-06 | 2.6926e-11 |
| 1 | 1.0131318160762 | 1.0131318154099 | 1.0197916823397 | 1.0197916823397 | 5.074e-06 | 2.6926e-11 |

Maximum horizontal error: 7.4162e-05 m.
Vertical position remains exactly 0 m. The particle remained wet, with constant
density, radius and condition, no degradation, and zero stochastic velocity.

All NetCDF values match the prescribed currents (0.2,0.3,0), wind (10,15),
and Stokes drift (0.2,0.4) m/s. Applying coefficients 0.02 and 0.03 gives
total horizontal velocity (0.406,0.612) m/s. The recorded current is checked
after the first interpolation; the age-zero VTU contains initialized zeros.

The closed-form longitude was independently checked against 32-point Gaussian
quadrature within 2e-15 degrees at both checkpoints. The longitude error is
consistent with forward Euler evaluating the changing spherical metric at
the beginning of each step. The latitude discrepancy is consistent with roundoff.

This case verifies simultaneous additive advection, windage and prescribed
Stokes drift. The paper buoyancy kernel is active at zero density contrast.
Nonzero settling/rising, stochastic transport and boundary interactions are
covered separately. The kernel sums transport contributions before integration;
this run does not establish the accuracy of a sequential operator-splitting method.

The generated LaTeX table, figure, CSV and JSON summary are in
`CoupledProcesses_Def_out`. Running the shell script replaces that folder.
