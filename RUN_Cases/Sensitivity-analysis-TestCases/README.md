# Arousa 3D sensitivity analysis

Forty-five Lagrangian cases on the same Cambados A1 source, the same
3-month window (`2024-02-23` to `2024-05-23`) and the same baseline as
`Arousa_3D_Def.xml` (original paper S8 / 1 mm, 1050 kg m⁻³):

- `DiffusionCoeffHorz = 1.0`, `DiffusionCoeffVert = 1.0e-5`
- `WindDragCoeff = 0`
- no wind / wave files
- degradation off (`Lin0 = Exp0 = Exp1 = Exp2 = 0`)
- resuspension off (`ResuspensionProb = 0`)

Only the axis named in the file changes. `OutputWriteTime` equals the 90-day
run length, so each case writes the start and the last instant used by the
paper tables.

## Run

With the model compiled and Conda available, from this folder:

```bash
./RunCase_S08-MaterialType.sh
```

Scripts locate the repository from their own path, so they can also be
launched from another directory. Outputs remain inside `${name}_out`.
Each script clears that folder when run again. Shared inputs are
`materialTypes.xml`, `ncNamesLibrary.xml`, `data/` and `nc_fields/`
(hydrodynamics and water properties from 22 Feb to 23 May 2024).

To regenerate every XML and script from this file:

```bash
python generate_cases.py
```

That does not touch existing `_out` folders.

## Cases

| ID | File | What changes |
|---:|---|---|
| 01 | S01-MaterialType | radius 0.2 mm, density 950 kg m-3 |
| 02 | S02-MaterialType | radius 0.2 mm, density 1050 kg m-3 |
| 03 | S03-MaterialType | radius 0.2 mm, density 1070 kg m-3 |
| 04 | S04-MaterialType | radius 0.6 mm, density 950 kg m-3 |
| 05 | S05-MaterialType | radius 0.6 mm, density 1050 kg m-3 |
| 06 | S06-MaterialType | radius 0.6 mm, density 1070 kg m-3 |
| 07 | S07-MaterialType | radius 1 mm, density 950 kg m-3 |
| 08 | S08-MaterialType | radius 1 mm, density 1050 kg m-3 |
| 09 | S09-MaterialType | radius 1 mm, density 1070 kg m-3 |
| 10 | S10-MaterialType | radius 2 mm, density 950 kg m-3 |
| 11 | S11-MaterialType | radius 2 mm, density 1050 kg m-3 |
| 12 | S12-MaterialType | radius 2 mm, density 1070 kg m-3 |
| 13 | S13-Degradation | mu0=0.14, mu2=0.040 1/day (Exp1=0.514) |
| 14 | S14-Degradation | mu0=0.14, mu2=0.40 1/day (Exp1=0.514) |
| 15 | S15-ResuspensionMethod1 | radius 0.6 mm, density 1050 kg m-3, tau_cr1=0.01, tau_cr2=0.10, eps0=0.0025 m |
| 16 | S16-ResuspensionMethod1 | radius 1 mm, density 1050 kg m-3, tau_cr1=0.01, tau_cr2=0.10, eps0=0.0025 m |
| 17 | S17-ResuspensionMethod1 | radius 2 mm, density 1050 kg m-3, tau_cr1=0.01, tau_cr2=0.10, eps0=0.0025 m |
| 18 | S18-ResuspensionMethod1 | radius 0.6 mm, density 1050 kg m-3, tau_cr1=0.01, tau_cr2=0.10, eps0=0.0050 m |
| 19 | S19-ResuspensionMethod1 | radius 1 mm, density 1050 kg m-3, tau_cr1=0.01, tau_cr2=0.10, eps0=0.0050 m |
| 20 | S20-ResuspensionMethod1 | radius 2 mm, density 1050 kg m-3, tau_cr1=0.01, tau_cr2=0.10, eps0=0.0050 m |
| 21 | S21-ResuspensionMethod1 | radius 0.6 mm, density 1050 kg m-3, tau_cr1=0.01, tau_cr2=0.10, eps0=0.0100 m |
| 22 | S22-ResuspensionMethod1 | radius 1 mm, density 1050 kg m-3, tau_cr1=0.01, tau_cr2=0.10, eps0=0.0100 m |
| 23 | S23-ResuspensionMethod1 | radius 2 mm, density 1050 kg m-3, tau_cr1=0.01, tau_cr2=0.10, eps0=0.0100 m |
| 24 | S24-ResuspensionMethod1 | radius 0.6 mm, density 1050 kg m-3, tau_cr1=0.01, tau_cr2=0.04, eps0=0.0025 m |
| 25 | S25-ResuspensionMethod1 | radius 1 mm, density 1050 kg m-3, tau_cr1=0.01, tau_cr2=0.04, eps0=0.0025 m |
| 26 | S26-ResuspensionMethod1 | radius 2 mm, density 1050 kg m-3, tau_cr1=0.01, tau_cr2=0.04, eps0=0.0025 m |
| 27 | S27-ResuspensionMethod1 | radius 0.6 mm, density 1050 kg m-3, tau_cr1=0.01, tau_cr2=0.04, eps0=0.0050 m |
| 28 | S28-ResuspensionMethod1 | radius 1 mm, density 1050 kg m-3, tau_cr1=0.01, tau_cr2=0.04, eps0=0.0050 m |
| 29 | S29-ResuspensionMethod1 | radius 2 mm, density 1050 kg m-3, tau_cr1=0.01, tau_cr2=0.04, eps0=0.0050 m |
| 30 | S30-ResuspensionMethod1 | radius 0.6 mm, density 1050 kg m-3, tau_cr1=0.01, tau_cr2=0.04, eps0=0.0100 m |
| 31 | S31-ResuspensionMethod1 | radius 1 mm, density 1050 kg m-3, tau_cr1=0.01, tau_cr2=0.04, eps0=0.0100 m |
| 32 | S32-ResuspensionMethod1 | radius 2 mm, density 1050 kg m-3, tau_cr1=0.01, tau_cr2=0.04, eps0=0.0100 m |
| 33 | S33-ResuspensionMethod2 | radius 0.6 mm, density 1050 kg m-3, variable tau_cr (Shields/Soulsby), eps0=0.0025 m |
| 34 | S34-ResuspensionMethod2 | radius 1 mm, density 1050 kg m-3, variable tau_cr (Shields/Soulsby), eps0=0.0025 m |
| 35 | S35-ResuspensionMethod2 | radius 2 mm, density 1050 kg m-3, variable tau_cr (Shields/Soulsby), eps0=0.0025 m |
| 36 | S36-ResuspensionMethod2 | radius 0.6 mm, density 1050 kg m-3, variable tau_cr (Shields/Soulsby), eps0=0.0050 m |
| 37 | S37-ResuspensionMethod2 | radius 1 mm, density 1050 kg m-3, variable tau_cr (Shields/Soulsby), eps0=0.0050 m |
| 38 | S38-ResuspensionMethod2 | radius 2 mm, density 1050 kg m-3, variable tau_cr (Shields/Soulsby), eps0=0.0050 m |
| 39 | S39-ResuspensionMethod2 | radius 0.6 mm, density 1050 kg m-3, variable tau_cr (Shields/Soulsby), eps0=0.0100 m |
| 40 | S40-ResuspensionMethod2 | radius 1 mm, density 1050 kg m-3, variable tau_cr (Shields/Soulsby), eps0=0.0100 m |
| 41 | S41-ResuspensionMethod2 | radius 2 mm, density 1050 kg m-3, variable tau_cr (Shields/Soulsby), eps0=0.0100 m |
| 42 | S42-Diffusion | Dh=1 m2/s, Dv=1e-4 m2/s |
| 43 | S43-Diffusion | Dh=1 m2/s, Dv=1e-6 m2/s |
| 44 | S44-Diffusion | Dh=2 m2/s, Dv=1e-5 m2/s |
| 45 | S45-Diffusion | Dh=4 m2/s, Dv=1e-5 m2/s |
