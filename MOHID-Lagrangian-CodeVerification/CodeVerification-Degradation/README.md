# Detritus degradation

One detritus tracer is followed for one hour at constant temperature 25 degrees C.
The current XML uses nonzero degradation coefficients and L(25 C)=1 to double
precision. Mass, volume and radius decrease with time. This removes the initial
temperature-limiter discontinuity while retaining a nonconstant reference.
The LaTeX table uses an Initial row at 0 h, then Analytical and Numerical
rows at 0.5 and 1 h. Columns contain mass, volume, radius and their relative
percentage errors. Euler and the
output interval both use 60 seconds.

All transport contributions are zero:

- The hydrodynamic NetCDF has u = v = w = 0 everywhere.
- No wind or wave inputs are supplied; their coefficients are zero.
- All diffusion coefficients and resuspension probability are zero.
- Particle and constant water densities both equal 1026 kg/m3, giving zero buoyancy.

Particle age still increases because the degradation law depends on time.
The tracer starts at
(1 degree, 1 degree, -100 m), inside an all-wet domain from 0.9 to 1.1 degrees
and -200 to 0 m. Temperature and salinity are initialized in materialTypes.xml;
the analyzer checks they remain 25 and 36 throughout the run.

Initial radius is 1 mm. Initial volume is 4*pi*r0^3/3 and initial mass is rho*V0.
Density is constant. The minimum-radius mass floor is included in the analytical
reference but is not reached during this test.

## Run

The input is included. Regenerate it when needed:

```bash
python analyze_results.py --prepare-inputs
```

Run and analyze separately in the MOHID-Lagrangian environment:

```bash
bash RunCase_Degradation.sh
python analyze_results.py
```

The shell script recreates `Degradation_Def_out`. Analysis never runs the model.
It writes `verification_table.tex`, `Degradation_verification.pdf`/`.png`,
`degradation.csv` and `summary.json` inside that folder, even when a numerical
check fails. `--input /path/to/output` analyzes another completed run.

This case checks the Euler implementation against the analytical law.
It does not establish compatibility with Euler2 or RK4.

## Reference and units

The temperature limiter follows the supplied BK1/BK2/BK3/BK4 expressions.
The current settings use ToptBMin=15, ToptBMax=30, TBMin=5, TBMax=35,
BK1=BK4=0.5 and BK2=BK3=0.999999999999. These give L(25 C)=1 to
double precision. The limiter remains temperature-dependent; it is not
identically one for all temperatures.
XML rate inputs are:

| Input | XML value | Units |
|---|---:|---|
| Lin0_DegradationRate | 0.03 | 1/day |
| Exp0_DegradationRate | 0.14 | dimensionless mixture fraction |
| Exp1_DegradationRate | 0.514 | 1/day |
| Exp2_DegradationRate | 0.04 | 1/day |

The model's setters divide the three rates by 86400. The Python reference
makes the same unit conversion before evaluating time in seconds:

```text
F(t) = (1 - L0*t) * [E0*exp(-E1*t) + (1-E0)*exp(-E2*t)]
M(t) = max(M0 * L(T) * F(t), rho * 4*pi*rmin^3/3)  for t > 0
V(t) = M(t)/rho
r(t) = [3*V(t)/(4*pi)]^(1/3)
```

At age zero the constructor initializes M0, V0 and r0 directly. The multiplier
L(T) is first applied during the first degradation update. Thus the current
formula has an initial reduction when L(T) differs from 1; the reference
preserves this behavior rather than normalizing the limiter away.

There is no independently stored mass field in this tracer. Numerical mass is
reconstructed as rho times the numerical VTU volume. Radius is compared as a
radius, not as diameter; diameter would be 2*r.

## Radius and volume

The degradation kernel returns rates for volume and radius. For Euler,
these rates advance both properties to their prescribed end-of-step values.
The reference compares mass, volume and radius at the same particle age.

Acceptance requires constant density and initial volume, one wet stationary
particle, no nonzero resolved or stochastic velocity, exact output checkpoints,
and mass, volume and same-time radius relative errors <= 1e-10.
Variable-temperature forcing and other integrators require separate checks.
