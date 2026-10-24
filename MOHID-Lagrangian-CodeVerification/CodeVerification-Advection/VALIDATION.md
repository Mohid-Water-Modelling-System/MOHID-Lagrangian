# Validation

Euler2 was rerun on 23 October 2026 after recompiling with
`mstime = time + dt` in `runStepMSEuler`. All four full-model runs passed
the particle identity, water-state, boundary and 3600 s age checks. Euler1
and RK4 rows below retain their previous full-model results.

| Method | Time step (s) | Position error (m) | Observed order |
|---|---:|---:|---:|
| Eluer1 | 120 | 2.67432 | 0.985 |
| Eluer1 | 60 | 1.35156 | 0.992 |
| Eluer1 | 30 | 0.679461 | 0.996 |
| Eluer1 | 15 | 0.340661 | — |
| Eluer2 | 120 | 0.00781116 | 2.010 |
| Eluer2 | 60 | 0.00193885 | 2.005 |
| Eluer2 | 30 | 0.000482909 | 2.003 |
| Eluer2 | 15 | 0.000120498 | — |
| RK4 | 120 | 3.34069e-07 | 3.997 |
| RK4 | 60 | 2.09295e-08 | 3.991 |
| RK4 | 30 | 1.31613e-09 | 3.860 |
| RK4 | 15 | 9.06239e-11 | — |

Euler1 approaches order one, corrected Euler2 approaches order two, and RK4
approaches order four. The finest RK4 pair is affected slightly by
interpolation and floating-point error. These are deterministic advection
results; they do not establish stochastic convergence order.

The Euler2 errors agree with an independent Heun calculation of the
analytical Cartesian field to within 1e-9 m for all four time steps.
The previous half-step corrector time gave orders 1.004, 1.002 and 1.001.
Its full-model results and the new executable SHA-256 are retained in
`build/investigations/advection_before_heun_retest_20261023_190620/`.

The final analytical position is approximately `(165.53293632, 128.65960373,
-46.95408390)` m. No particle reached the domain boundaries.

## Grid check

RK4 was also run with 401 latitude points instead of 201, keeping the same
longitude/depth/time coordinates. The values below record that earlier grid check. Its investigation output
folder is no longer present after the rebuild; the grid check was not rerun
for this Euler2 correction.

| Time step (s) | Original error (m) | Refined error (m) | Change in endpoint (m) |
|---|---:|---:|---:|
| 120 | 3.34069e-07 | 3.34063e-07 | 1.32729e-11 |
| 15 | 9.06239e-11 | 8.38e-11 | 1.30456e-11 |

The finest-step endpoint changes by about 13 picometres, but this is
a noticeable fraction of its exceptionally small error. Interpret the finest
observed order with that interpolation floor in mind.

The 12 XML files parse and the 12 shell scripts pass `bash -n`. The shell
Conda activation was not exercised here; the same preprocessor/model commands
were run directly with the project Python environment.

## Configuration checks

All 12 source configurations pass the shared validation used by the run
scripts and analysis. A deliberately incorrect Start time was rejected.
Running its shell script left a pre-existing output file intact. The
`--sync-cases` command restored the correct settings while preserving the
XML comments. Existing VTU analysis also completes with these checks enabled.
