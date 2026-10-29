# Turbulence tests

Prepared by Mohsen Shabani | Email: shabani.mohsen@outlook.com

Two cases with 20,000 passive `base` particles in a 20 × 20 km water domain,
200 m deep. Particles start at the center, 100 m below the surface. The NetCDF
contains zero current in all directions. No wind or wave files are used.

## Run

From this folder, with the model already compiled:

```bash
./RunCase_RandomWalk.sh
./RunCase_RandomFlight.sh
```

Each script creates its own output folder:

- `RandomWalk_Def_out`
- `RandomFlight_Def_out`

Scripts locate the repository from their own path and can be launched from
any directory. Outputs remain inside the case folder.
The scripts clear the corresponding output folder before each run, as in the
other Mohsen cases. Keep a copy of any results you want to retain.
Simulation settings are in the corresponding XML file. The scripts activate the `MOHID-Lagrangian` Conda environment.
The `data` folder is copied from NewParam, with the output fields needed for
these statistical tests in `data/outputFields.xml`.

## Analyze

After the runs finish:

```bash
conda activate MOHID-Lagrangian
python analyze_results.py
```

Use `--case rw` or `--case rf` to analyze one case. Tables, figures, statistics
and the pass/fail report are saved in the same `_out` folder:

- `statistics.csv`
- `verification_table.tex` (reference, numerical value and relative error in percent)
- `*_MSD_verification.png` and `.pdf`
- `summary.json`

If the NetCDF input is missing, recreate it before running:

```bash
python analyze_results.py --prepare-inputs
```

To redraw only the figures from the existing statistics:

```bash
python analyze_results.py --figures-only
```

Figures use serif fonts, open markers and a black dashed analytical reference.
PNG files are saved at 600 dpi, with vector PDF copies.

Relative error is `100 * abs(numerical - reference) / abs(reference)`.
The LaTeX tables follow the manuscript rows (six for Random Walk and eight
for Random Flight). For a zero reference, relative error is `--` and the
numerical estimate includes one standard error. The tables require `booktabs`
and `amsmath`.

Python packages are listed in `requirements.txt`.

## Settings

Both tests use Euler integration and a physical duration of 3600 s.

| Setting | Random Walk | Random Flight |
|---|---:|---:|
| Time step | 10 s | 5 s |
| Horizontal diffusivity | 2 m²/s | — |
| Vertical diffusivity | 1e-5 m²/s | — |
| RMS turbulent speed | — | 0.1 m/s |
| Mixing length | — | 6 m |
| Eddy lifetime | — | 60 s |

The XML end time includes two extra time steps for the solver's end condition.
The analysis uses particle `age` and checks 3600 s of motion. The source emits
only once; its initial emission counter means `19999` gives 20,000 particles.

The analysis compares VTU statistics with the analytical mean, variance and MSD.
Random Flight also checks velocity persistence and renewal timing. Acceptance
limits use six estimated standard errors. Numerical results are in
[VALIDATION.md](VALIDATION.md).
