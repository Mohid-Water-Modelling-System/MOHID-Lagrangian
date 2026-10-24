# MOHID-Lagrangian code verification

Fundamental verification cases for transport and water properties.

- [Buoyancy-driven motion](CodeVerification-Buoyancy-deriven/README.md): constant-property
  settling and rising trajectories at 0.5 and 1 hour.
- [Detritus degradation](CodeVerification-Degradation/README.md): mass, volume
  and radius at constant temperature, with all transport contributions zero.
- [Coupled processes](CodeVerification-CoupledProcesses-Adv+Wd+Wv/README.md):
  simultaneous advection, windage and prescribed Stokes drift.
- [Stokes drift](CodeVerification-StokesDrift/README.md): prescribed surface
  drift, coefficient scaling and spherical coordinate conversion.
- [Wind-induced drift](CodeVerification-WindDrift/README.md): constant-wind
  surface transport and spherical coordinate conversion.
- [Advection](CodeVerification-Advection/README.md): analytical trajectories
  and temporal convergence of Euler, second-order Euler and RK4.
- [Water properties](CodeVerification-WaterProperties/README.md): independent
  density and dynamic-viscosity selections, checked against the model library.
- [Turbulence](CodeVerification-Turbulence/README.md): random-walk and
  random-flight statistical checks under zero resolved current.

Shell scripts resolve their own case folder and can be launched from any
directory. Outputs are saved beside the case XML. Analyze
existing outputs with its `analyze_results.py`; analysis does not run the
model. Each case keeps its input, configuration and output files together.
