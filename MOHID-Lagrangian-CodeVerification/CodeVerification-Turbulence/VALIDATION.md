# Test results

Prepared by Mohsen Shabani | Email: shabani.mohsen@outlook.com

Run on 23 October 2026 with 20,000 particles, one thread and one block.
Both tests cover 3600 s of particle motion.

| Quantity | Analytical | Measured |
|---|---:|---:|
| Random Walk mean x displacement (m) | 0 | 0.856430 |
| Random Walk mean y displacement (m) | 0 | -0.168558 |
| Random Walk horizontal MSD (m²) | 28800 | 29070.79 |
| Random Walk x variance (m²) | 14400 | 14441.13 |
| Random Walk y variance (m²) | 14400 | 14630.35 |
| Random Walk z variance (m²) | 0.072 | 0.072605 |
| Random Walk fitted diffusivity (m²/s) | 2 | 2.02968 |
| Random Flight mean x velocity (m/s) | 0 | -0.00107245 |
| Random Flight mean y velocity (m/s) | 0 | -0.00064103 |
| Random Flight horizontal MSD (m²) | 2160 | 2159.58 |
| Random Flight x velocity variance (m²/s²) | 0.005 | 0.004977 |
| Random Flight y velocity variance (m²/s²) | 0.005 | 0.004973 |
| Random Flight velocity covariance (m²/s²) | 0 | 5.8791775e-06 |
| Random Flight RMS speed (m/s) | 0.1 | 0.099754 |
| Random Flight renewal interval (s) | 60 | 60 |
| Random Flight fitted diffusivity (m²/s) | 0.15 | 0.150924 |

Both tests passed the statistical checks. All 20,000 particles remained in
water with zero resolved current. The minimum boundary clearance was 98.81 m
for Random Walk and 100 m for Random Flight.

Random Flight retained its velocity between renewals and renewed every 60 s.
Its vertical displacement stayed zero. MSD differences were within 2.3 standard
errors for both tests; the acceptance limit was six standard errors.

VTUs, tables and figures are in the corresponding `_out` folders. These values
refer to this run; another run may differ because of random sampling.
