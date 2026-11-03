#!/usr/bin/env python3
"""Generate the S01-S45 XML + RunCase.sh suite for the Arousa sensitivity analysis.

All cases share the same 3-month window (2024-02-23 to 2024-05-23), the same
source (CambadosA1) and the same baseline physics as the S08 template. Only
the axis under test changes. Naming follows CodeVerification-Advection:
Sxx-AnalysisType.xml and RunCase_Sxx-AnalysisType.sh.
"""
from __future__ import annotations

from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent

START = "2024 02 23 00 00 00"
END = "2024 05 23 00 00 00"
OWT = "7776000"  # 90 days: start + end only (tables use the last instant)

# Kelp_detritus_i_j : i = radius class (1..4), j = density class (1..3)
MATERIALS = {
    "1_1": {"density": 950.0, "radius": 0.0002, "volume": "0.0000000000335", "area": "0.00000050265"},
    "1_2": {"density": 1050.0, "radius": 0.0002, "volume": "0.0000000000335", "area": "0.00000050265"},
    "1_3": {"density": 1070.0, "radius": 0.0002, "volume": "0.0000000000335", "area": "0.00000050265"},
    "2_1": {"density": 950.0, "radius": 0.0006, "volume": "0.0000000009052", "area": "0.0000045239"},
    "2_2": {"density": 1050.0, "radius": 0.0006, "volume": "0.0000000009052", "area": "0.0000045239"},
    "2_3": {"density": 1070.0, "radius": 0.0006, "volume": "0.0000000009052", "area": "0.0000045239"},
    "3_1": {"density": 950.0, "radius": 0.001, "volume": "0.000000004189", "area": "0.00001257"},
    "3_2": {"density": 1050.0, "radius": 0.001, "volume": "0.000000004189", "area": "0.00001257"},
    "3_3": {"density": 1070.0, "radius": 0.001, "volume": "0.000000004189", "area": "0.00001257"},
    "4_1": {"density": 950.0, "radius": 0.002, "volume": "0.00000003351", "area": "0.00005027"},
    "4_2": {"density": 1050.0, "radius": 0.002, "volume": "0.00000003351", "area": "0.00005027"},
    "4_3": {"density": 1070.0, "radius": 0.002, "volume": "0.00000003351", "area": "0.00005027"},
}

RADIUS_MM = {"1": 0.2, "2": 0.6, "3": 1.0, "4": 2.0}
DENSITY = {"1": 950, "2": 1050, "3": 1070}


def default_physics() -> dict:
    return {
        "property": "Kelp_detritus_3_2",
        "dh": "1.0",
        "dv": "1.0e-5",
        "resusp_method": "1",
        "resusp_prob": "0",
        "tau1": "100",
        "tau2": "100",
        "rugosity": "0.0025",
        "lin0": "0.00000",
        "exp0": "0.0",
        "exp1": "0.0",
        "exp2": "0.0",
    }


def material_type_cases() -> list[dict]:
    cases = []
    sid = 1
    for i in ("1", "2", "3", "4"):
        for j in ("1", "2", "3"):
            p = default_physics()
            p["property"] = f"Kelp_detritus_{i}_{j}"
            cases.append(
                {
                    "id": sid,
                    "kind": "MaterialType",
                    "title": (
                        f"radius {RADIUS_MM[i]:g} mm, "
                        f"density {DENSITY[j]} kg m-3"
                    ),
                    **p,
                }
            )
            sid += 1
    return cases


def degradation_cases() -> list[dict]:
    c13 = default_physics()
    c13.update(
        id=13,
        kind="Degradation",
        title="mu0=0.14, mu2=0.040 1/day (Exp1=0.514)",
        exp0="0.14",
        exp1="0.514",
        exp2="0.040",
    )
    c14 = default_physics()
    c14.update(
        id=14,
        kind="Degradation",
        title="mu0=0.14, mu2=0.40 1/day (Exp1=0.514)",
        exp0="0.14",
        exp1="0.514",
        exp2="0.40",
    )
    return [c13, c14]


def resuspension_method1_cases() -> list[dict]:
    radii = (("2_2", 0.6), ("3_2", 1.0), ("4_2", 2.0))
    cases = []
    sid = 15
    for tau2 in ("0.10", "0.04"):
        for eps0 in ("0.0025", "0.0050", "0.0100"):
            for prop, rmm in radii:
                p = default_physics()
                p.update(
                    id=sid,
                    kind="ResuspensionMethod1",
                    title=(
                        f"radius {rmm:g} mm, density 1050 kg m-3, "
                        f"tau_cr1=0.01, tau_cr2={tau2}, eps0={eps0} m"
                    ),
                    property=f"Kelp_detritus_{prop}",
                    resusp_method="1",
                    resusp_prob="100",
                    tau1="0.01",
                    tau2=tau2,
                    rugosity=eps0,
                )
                cases.append(p)
                sid += 1
    return cases


def resuspension_method2_cases() -> list[dict]:
    radii = (("2_2", 0.6), ("3_2", 1.0), ("4_2", 2.0))
    cases = []
    sid = 33
    for eps0 in ("0.0025", "0.0050", "0.0100"):
        for prop, rmm in radii:
            p = default_physics()
            p.update(
                id=sid,
                kind="ResuspensionMethod2",
                title=(
                    f"radius {rmm:g} mm, density 1050 kg m-3, "
                    f"variable tau_cr (Shields/Soulsby), eps0={eps0} m"
                ),
                property=f"Kelp_detritus_{prop}",
                resusp_method="2",
                resusp_prob="100",
                tau1="0.01",
                tau2="0.04",
                rugosity=eps0,
            )
            cases.append(p)
            sid += 1
    return cases


def diffusion_cases() -> list[dict]:
    specs = (
        (42, "Dh=1 m2/s, Dv=1e-4 m2/s", {"dv": "1.0e-4"}),
        (43, "Dh=1 m2/s, Dv=1e-6 m2/s", {"dv": "1.0e-6"}),
        (44, "Dh=2 m2/s, Dv=1e-5 m2/s", {"dh": "2.0"}),
        (45, "Dh=4 m2/s, Dv=1e-5 m2/s", {"dh": "4.0"}),
    )
    cases = []
    for sid, title, extra in specs:
        p = default_physics()
        p.update(id=sid, kind="Diffusion", title=title, **extra)
        cases.append(p)
    return cases


def all_cases() -> list[dict]:
    cases = (
        material_type_cases()
        + degradation_cases()
        + resuspension_method1_cases()
        + resuspension_method2_cases()
        + diffusion_cases()
    )
    assert [c["id"] for c in cases] == list(range(1, 46))
    return cases


def case_name(case: dict) -> str:
    return f"S{case['id']:02d}-{case['kind']}"


def write_xml(case: dict) -> None:
    name = case_name(case)
    comment = (
        f"Sensitivity analysis {name}: {case['title']}. "
        "Shared window 2024-02-23 to 2024-05-23 (90 days). "
        "Template is Arousa_3D_Def.xml (original S8 / M32)."
    )
    xml = f"""<?xml version="1.0" encoding="UTF-8" ?>
<case>
	<!-- {comment} -->
	<execution>
		<parameters>
			<parameter key="Start" value="{START}" comment="Date of initial instant" units_comment="space delimited ISO 8601 format up to seconds" />
			<parameter key="End"   value="{END}" comment="Date of final instant" units_comment="ISO format" />
			<parameter key="BaseDateTime" value="1950 01 01 00 00 00" comment="Base Date for time stamping purposes. Optional, default = 1950/1/1" units_comment="ISO format" />
			<parameter key="Integrator" value="1" comment="Integration Algorithm 1:Euler, 2:Multi-Step Euler, 3:RK4 (default=1)" />
			<parameter key="Threads" value="5" comment="Computation threads for shared memory computation (default=auto)" />
			<parameter key="BufferSize" value="86400" comment="Optional parameter. Controls input frequency" units_comment="seconds" />
			<parameter key="OutputWriteTime" value="{OWT}" comment="Write interval = 90-day run length. Start + end only (paper tables use the last instant)." units_comment="seconds" />
			<parameter key="OutputFormat" value="2" comment="Output file format. NetCDF=1; VTK=2 (default=2)" />
		</parameters>
		<variableNaming>
			<file name="ncNamesLibrary.xml"/>
		</variableNaming>
		<outputFields>
			<file name="data/outputFields.xml"/>
		</outputFields>
	</execution>
	<caseDefinitions>
		<inputData>
			<inputDataDir name="nc_fields/hydrodynamic" type="hydrodynamic"/>
      		<inputDataDir name="nc_fields/waterProperties" type="waterProperties"/>
		</inputData>
		<simulation>
			<resolution dp="300" units_comment="metres (m)"/>
			<timestep dt="60" units_comment="seconds (s)"/>
			<BoundingBoxMin x="-9.10" y="42.39" z="-140" units_comment="(deg,deg,m)"/>
			<BoundingBoxMax x="-8.72" y="42.68" z="3" units_comment="(deg,deg,m)"/>
            <VerticalVelMethod value="1" comment="1:From velocity fields, 2:Divergence based, 3:Disabled. Default = 1" />
			<BathyminNetcdf value="1" comment="bathymetry is a property in the netcdf. 1:true, 0:false (computes from layer depth and openPoints. Default = 1"/>
			<RugosityinNetcdf value="0" comment="Rugosity of the seabed is a property in the netcdf. 1:true, 0:false (Default = 0, read from constant input value)"/>
			<D50inNetcdf value="0" comment="D50 (median diameter) of the seabed is a property in the netcdf. 1:true, 0:false (Default = 0, read from constant input value)"/>
			<RemoveLandTracer value="0" comment="Remove tracers on land 0:No, 1:Yes. Default = 0" />
            <ResuspensionCriticalShearMethod value="{case['resusp_method']}" comment="1:Constant tau_cr, 2:Shields/Soulsby variable tau_cr." />
        </simulation>
		<sourceDefinitions>
			<source>
				<setsource id="1" name="CambadosA1" />
				<resolution x="20" y="20" z="150" units_comment="metres (m)"/>
				<rate_seconds value="900" comment="emission step in seconds. 3600 is a tracer per hour. it should be bigger than dt_simulation" />
                <rate_trcPerEmission value="1" comment="number of tracers emiited every rate_seconds. 5 is a 5 tracers per rate_seconds" />
				<active start="0" end="1800" comment="example: start='12.7' end='end'; start='0.0' end='95' " units_comment="seconds (s)" />
				<polygon>
					<file name="data/CambadosA1.kmz"/>
                    <verticalBoundingBox min="-1.0" max="0.0"/>
				</polygon>
			</source>
		</sourceDefinitions>
		<sourceTypes>
			<types>
                <type source="1" type="detritus" property="{case['property']}" comment="" />
			</types>
			<file name="materialTypes.xml"/>
		</sourceTypes>
		<constants>
			<DiffusionCoeffHorz value="{case['dh']}" comment="Horizontal diffusion coefficient. Default = 1.0" units_comment="m2/s" />
			<DiffusionCoeffVert value="{case['dv']}" comment="Vertical diffusion coefficient; used directly. Default = 0.000001." units_comment="m2/s" />
			<Rugosity value="{case['rugosity']}" comment="Rugosity / bed roughness length eps0. Default = 0.0025" units_comment="m" />
			<D50 value="0.0025" comment="D50 (Median diameter of seabed) value. Default = 0.0025" units_comment="m" />
			<D50Density value="2800.0" comment="D50 density (Median diameter of seabed density) value. Default = 2500.0" units_comment="kg/m3" />
			<BeachingLevel value="-3.0" comment="Level above which beaching can occur. Default = -3.0" units_comment="m" />
			<BeachingStopProb value="100" comment="Probablity of beaching stopping a tracer. Default = 50%" units_comment="%" />
			<WindDragCoeff value="0.0" comment="Windage coefficient. Default = 0.0" units_comment="" />

			<BedLoadThickness value="0.04" comment="BedLoad thickness. Default = 0.04" units_comment="m" />
			<ResuspensionProb value="{case['resusp_prob']}" comment="Probablity of Resuspension of a tracer. Default = 0%" units_comment="%" />
			<ResuspensionResidenceTime value="8640000.0" comment="Resuspension Residence Time. Default = 10*24*60*60" units_comment="sec" />
			<ResuspensionCriticalShear1 value="{case['tau1']}" comment="Resuspension Critical Shear 1. Unused if method=2 or ResuspensionProb=0." units_comment="N/m2" />
			<ResuspensionCriticalShear2 value="{case['tau2']}" comment="Resuspension Critical Shear 2. Unused if method=2 or ResuspensionProb=0." units_comment="N/m2" />

			<ToptBMin value="24.8" comment="minimum temperature of the optimal interval for the Bacteria growth. Default = 24.8" units_comment="degrees" />
            <ToptBMax value="25.1" comment="maximum temperature of the optimal interval for the Bacteria growth. Default = 25.1" units_comment="degrees" />
            <TBMin value="5" comment="minimum tolerable temperature of the interval for the Bacteria growth. Default = 5" units_comment="degrees" />
            <TBMax value="35" comment="maximum tolerable temperature of the interval for the Bacteria growth. Default = 35" units_comment="degrees" />
            <BK1 value="0.9999" comment="constant to control temperature response curve shape. Default = 0.05" units_comment="%" />
            <BK2 value="0.9999" comment="constant to control temperature response curve shape. Default = 0.98" units_comment="%" />
            <BK3 value="0.9999" comment="constant to control temperature response curve shape. Default = 0.98" units_comment="%" />
            <BK4 value="0.9999" comment="constant to control temperature response curve shape. Default = 0.05" units_comment="%" />
			<Lin0_DegradationRate value="{case['lin0']}" comment="Lin0_DegradationRate. Default = 0.030/86400" units_comment="Input(1/day).In the code, 1/day converts 1/sec" />
			<Exp0_DegradationRate value="{case['exp0']}" comment="Exp0_DegradationRate. Default = 0.14" units_comment=" unit=1" />
			<Exp1_DegradationRate value="{case['exp1']}" comment="Exp1_DegradationRate. Default = 0.514/86400" units_comment="Input(1/day).In the code, 1/day converts 1/sec" />
			<Exp2_DegradationRate value="{case['exp2']}" comment="Exp2_DegradationRate. Default = 0.040/86400" units_comment="Input(1/day).In the code, 1/day converts 1/sec" />
		</constants>
	</caseDefinitions>
</case>
"""
    (CASE_DIR / f"{name}.xml").write_text(xml, encoding="utf-8")


def write_sh(case: dict) -> None:
    name = case_name(case)
    sh = f"""#!/bin/bash

# Resolve paths from this script so it can be launched from any directory.
caseDir="$(cd -- "$(dirname -- "${{BASH_SOURCE[0]}}")" && pwd)" || exit 1
repoRoot="$(cd -- "$caseDir/../../.." && pwd)" || exit 1
cd -- "$caseDir" || exit 1

clear

# Activate conda environment from shell. Minimum conda version required 4.6
eval "$(conda shell.bash hook)"
conda activate MOHID-Lagrangian

# Increase stacksizes in Linux.
ulimit -s hard
export OMP_STACKSIZE=512M

# "name" and "dirout" are named according to the testcase
name={name}
dirout=${{name}}_out
data=data

# "executables" are renamed and called from their directory
tools="${{repoRoot}}/build/bin"
mohidlagrangian=${{tools}}/MOHIDLagrangian

preprocessorDir="${{repoRoot}}/src/MOHIDLagrangianPreProcessor"
mohidPreprocessor=${{preprocessorDir}}/MOHIDLagrangianPreProcessor.py

postProcessorDir="${{repoRoot}}/src/MOHIDLagrangianPostProcessor"
mohidPostprocessor=${{postProcessorDir}}/MOHIDLagrangianPostProcessor.py

# "dirout" is created to store results or it is cleaned if it already exists
if [ -e "$dirout" ]; then
  rm -f -r "$dirout"
fi
mkdir "$dirout"

cp "${{name}}.xml" "$dirout/"
cp -r "$data" "$dirout/"

# CODES are executed according the selected parameters of execution in this testcase

python "$mohidPreprocessor" -i "$dirout/${{name}}.xml" -o "$dirout"
errcode=$?

if [ $errcode -eq 0 ]; then
  "$mohidlagrangian" -i "$dirout/${{name}}.xml" -o "$dirout"
  errcode=$?
fi

# python -W ignore "$mohidPostprocessor" -i "${{name}}.xml" -o "$dirout"

if [ $errcode -eq 0 ]; then
  echo All done
else
  echo Execution aborted
fi
read -n1 -r -p "Press any key to continue..." key
echo
"""
    path = CASE_DIR / f"RunCase_{name}.sh"
    path.write_text(sh, encoding="utf-8")
    path.chmod(0o755)


def write_material_types() -> None:
    lines = [
        '<?xml version="1.0" encoding="UTF-8" ?>',
        "<materials>",
        "	<detritus>",
    ]
    for key, mat in MATERIALS.items():
        i, j = key.split("_")
        lines.extend(
            [
                f"		<!-- r{i} = {RADIUS_MM[i]:g} mm, rho{j} = {DENSITY[j]} kg/m3 -->",
                f"		<Kelp_detritus_{key}>",
                '			<particulate value="1" />',
                f'			<density value="{mat["density"]}" />',
                f'			<radius value="{mat["radius"]}" />',
                f'			<volume value="{mat["volume"]}" />',
                f'			<area value="{mat["area"]}" />',
                '			<property key="condition" value="0.95" />',
                '			<property key="radius_cr_min" value="0.0000007" />',
                '			<property key="radius_cr_max" value="0.50" />',
                f"		</Kelp_detritus_{key}>",
            ]
        )
    lines.extend(["	</detritus>", "</materials>", ""])
    (CASE_DIR / "materialTypes.xml").write_text("\n".join(lines), encoding="utf-8")


def write_readme(cases: list[dict]) -> None:
    rows = []
    for c in cases:
        rows.append(f"| {c['id']:02d} | {case_name(c)} | {c['title']} |")
    body = f"""# Arousa 3D sensitivity analysis

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
launched from another directory. Outputs remain inside `${{name}}_out`.
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
{chr(10).join(rows)}
"""
    (CASE_DIR / "README.md").write_text(body, encoding="utf-8")


def main() -> None:
    cases = all_cases()
    write_material_types()
    for case in cases:
        write_xml(case)
        write_sh(case)
    write_readme(cases)
    print(f"Wrote {len(cases)} XML + {len(cases)} RunCase scripts in {CASE_DIR}")


if __name__ == "__main__":
    main()
