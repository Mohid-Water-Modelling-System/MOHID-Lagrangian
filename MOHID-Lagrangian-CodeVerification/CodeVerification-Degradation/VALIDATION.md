# Detritus degradation verification

The one-hour Euler run completed and passed all checks. One detritus particle
remained stationary while its mass, volume and radius decreased. The time step
and output interval are 60 s. The temperature is fixed at 25 C and L=1 to
double precision. All four degradation coefficients have the values in README.md.

The table compares the initial particle with analytical and numerical values
at particle ages 0.5 and 1 hour. One hour is sufficient to resolve the property
changes for these coefficients.

| Time (h) | Quantity | Analytical | Numerical | Relative error (%) |
|---|---|---:|---:|---:|
| 0.5 | Mass (kg) | 4.28553129072e-06 | 4.28553129072e-06 | 0 |
| 0.5 | Volume (m3) | 4.17693108257e-09 | 4.17693108257e-09 | 0 |
| 0.5 | Radius (m) | 0.000999055389031 | 0.000999055389031 | 0 |
| 1 | Mass (kg) | 4.27344642573e-06 | 4.27344642573e-06 | 0 |
| 1 | Volume (m3) | 4.16515246172e-09 | 4.16515246172e-09 | 0 |
| 1 | Radius (m) | 0.000998115418759 | 0.000998115418759 | 0 |

Maximum relative errors over the full trajectory (fractions, not percentages):

- Mass: 3.94814e-16
- Volume: 3.95585e-16
- Radius: 2.17243e-16

Maximum coordinate drift: 0. All resolved and stochastic velocities are zero.
The density and initial volume remain constant, and exactly one wet particle
is present throughout. The minimum-radius floor is not reached.

Mass is reconstructed from output density times volume. It is not an
independently stored state variable. Errors are calculated from unrounded data.

Outputs in Degradation_Def_out:

- verification_table.tex: initial/analytical/numerical particle table.
- Degradation_verification.pdf and .png: mass, volume, radius and relative errors.
- degradation.csv: full trajectory comparisons.
- summary.json: settings, numerical checks and checkpoint values.

RunCase_Degradation.sh runs preprocessing and the model. Run analyze_results.py
afterwards to generate the table and figures. Analysis never runs the model.

This verifies the specified Euler case at constant temperature. It does not
establish convergence order or accuracy under variable-temperature forcing.
