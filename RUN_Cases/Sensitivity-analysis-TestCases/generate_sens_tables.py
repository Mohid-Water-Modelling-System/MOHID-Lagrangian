#!/usr/bin/env python3
"""Build the sensitivity-analysis LaTeX tables from the final VTU of each case.

Classification (state is not used; the final output writes it as zero):

- beach (B): within 4 cm of the bed and local depth shallower than 3 m
- seabed (S): within 4 cm of the bed and local depth of 3 m or more
- ocean (O): more than 4 cm above the bed
- out of domain (OoD): emitted tracers absent from the final VTU

Zones follow the previous post-processing: Z1 is the Cambados A1 polygon,
Z2 is the 300 m buffer around it, and Z3 is everything else.

Degradation scenarios are written with a mass row and a number row.
Mass percentages use the emitted mass as the denominator.
"""
from __future__ import annotations

import re
import zipfile
from pathlib import Path

import numpy as np

CASE_DIR = Path(__file__).resolve().parent
OUT_TEX = CASE_DIR / "sensitivity_tables.tex"
OUT_TEX_COMPACT = CASE_DIR / "sensitivity_tables_2.tex"

BED_DISTANCE_M = 0.04
BEACH_DEPTH_M = 3.0
ZONE2_DISTANCE_M = 300.0
REFERENCE_ID = 8  # undegraded S08, included in the degradation table

NUMBER_SPEC = (
    r"l | *{3}{S[table-format=3.2]} S[table-format=3.2] "
    r"| *{3}{S[table-format=3.2]} S[table-format=3.2] "
    r"| *{3}{S[table-format=3.2]} S[table-format=3.2] "
    r"| *{3}{S[table-format=3.2]} S[table-format=3.2]"
)
MASS_SPEC = "l " + NUMBER_SPEC


def read_vtu(path: Path) -> dict[str, np.ndarray]:
    """Read one appended-data VTU. Coordinates come back as (N, 3)."""
    data = path.read_bytes()
    xml_end = data.find(b"<AppendedData")
    if xml_end < 0:
        raise ValueError(f"{path.name} has no appended data")
    xml = data[:xml_end].decode("latin1")
    payload = data[data.find(b"_", xml_end) + 1 :]
    piece_re = re.compile(
        r'<Piece NumberOfPoints="\+?(\d+)" NumberOfCells="\+?(\d+)">(.*?)</Piece>',
        re.S,
    )
    array_re = re.compile(
        r'<DataArray type="([^"]+)" NumberOfComponents="(\d+)" '
        r'Name="([^"]+)" format="appended" offset="(\d+)"/>'
    )
    dtypes = {
        "Float64": np.float64,
        "Float32": np.float32,
        "Int32": np.int32,
        "Int8": np.int8,
    }
    chunks: dict[str, list[np.ndarray]] = {}
    for piece in piece_re.finditer(xml):
        for match in array_re.finditer(piece.group(3)):
            dtype_name, ncomp, name, offset = match.groups()
            offset = int(offset)
            ncomp = int(ncomp)
            size = int.from_bytes(payload[offset : offset + 4], "little")
            raw = payload[offset + 4 : offset + 4 + size]
            arr = np.frombuffer(raw, dtype=dtypes[dtype_name])
            if ncomp > 1:
                arr = arr.reshape(-1, ncomp)
            chunks.setdefault(name, []).append(arr)
    if "Points" not in chunks:
        raise ValueError(f"{path.name} has no points")
    arrays = {name: np.concatenate(parts) for name, parts in chunks.items()}
    if "id" in arrays:
        _, unique = np.unique(arrays["id"], return_index=True)
        if len(unique) != len(arrays["id"]):
            arrays = {name: values[unique] for name, values in arrays.items()}
    return arrays


def tracer_files(folder: Path, stem: str) -> list[Path]:
    return sorted(folder.glob(f"{stem}_[0-9][0-9][0-9][0-9][0-9].vtu"))


def emission_count(xml_path: Path) -> int:
    text = xml_path.read_text()
    rate = float(re.search(r'rate_seconds value="([^"]+)"', text).group(1))
    active = re.search(r'<active start="([^"]+)" end="([^"]+)"', text)
    start = float(active.group(1))
    end = float(active.group(2))
    if rate <= 0:
        raise ValueError(f"rate_seconds in {xml_path.name} is {rate}")
    count = 0
    time = start
    while time < end - 1e-6:
        count += 1
        time += rate
    return count


def load_polygon(path: Path) -> np.ndarray:
    """Return unique ring vertices as (lon, lat). Accepts a KML saved as .kmz."""
    raw = path.read_bytes()
    if raw[:2] == b"PK":
        with zipfile.ZipFile(path) as archive:
            kml_name = next(name for name in archive.namelist() if name.endswith(".kml"))
            text = archive.read(kml_name).decode("utf-8")
    else:
        text = raw.decode("utf-8")
    coords = re.search(r"<coordinates>(.*?)</coordinates>", text, re.S).group(1)
    points = []
    for token in coords.split():
        lon, lat, *_ = token.split(",")
        points.append((float(lon), float(lat)))
    ring = np.asarray(points, dtype=float)
    if len(ring) > 1 and np.allclose(ring[0], ring[-1]):
        ring = ring[:-1]
    return ring


def to_local_xy(lon: np.ndarray, lat: np.ndarray, origin: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lat0 = np.deg2rad(origin[1])
    x = (lon - origin[0]) * (111_320.0 * np.cos(lat0))
    y = (lat - origin[1]) * 110_540.0
    return x, y


def points_in_polygon(x: np.ndarray, y: np.ndarray, poly: np.ndarray) -> np.ndarray:
    inside = np.zeros(x.shape, dtype=bool)
    px, py = poly[:, 0], poly[:, 1]
    j = len(poly) - 1
    for i in range(len(poly)):
        dy = py[j] - py[i]
        x_cross = (px[j] - px[i]) * (y - py[i]) / (dy if dy != 0 else 1.0) + px[i]
        crosses = ((py[i] > y) != (py[j] > y)) & (x < x_cross)
        inside ^= crosses
        j = i
    return inside


def distance_to_polygon(x: np.ndarray, y: np.ndarray, poly: np.ndarray) -> np.ndarray:
    best = np.full(x.shape, np.inf)
    for i in range(len(poly)):
        ax, ay = poly[i]
        bx, by = poly[(i + 1) % len(poly)]
        abx, aby = bx - ax, by - ay
        denom = abx * abx + aby * aby
        t = np.clip(((x - ax) * abx + (y - ay) * aby) / denom, 0.0, 1.0)
        best = np.minimum(best, np.hypot(x - (ax + t * abx), y - (ay + t * aby)))
    return best


def assign_zones(lon: np.ndarray, lat: np.ndarray, ring_lonlat: np.ndarray) -> np.ndarray:
    origin = ring_lonlat.mean(axis=0)
    poly = np.column_stack(to_local_xy(ring_lonlat[:, 0], ring_lonlat[:, 1], origin))
    x, y = to_local_xy(lon, lat, origin)
    inside = points_in_polygon(x, y, poly)
    zone = np.full(len(lon), 3, dtype=int)
    zone[inside] = 1
    outside = ~inside
    if np.any(outside):
        near = distance_to_polygon(x[outside], y[outside], poly) <= ZONE2_DISTANCE_M
        zone[np.flatnonzero(outside)[near]] = 2
    return zone


def particle_mass(radius: np.ndarray, density: np.ndarray) -> np.ndarray:
    return (4.0 / 3.0) * np.pi * radius**3 * density


def classify(z: np.ndarray, bathymetry: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    above_bed = z - bathymetry
    near_bed = above_bed < BED_DISTANCE_M
    shallow = (-bathymetry) < BEACH_DEPTH_M
    beach = near_bed & shallow
    seabed = near_bed & ~shallow
    ocean = ~near_bed
    return beach, seabed, ocean


def empty_block() -> dict[str, float]:
    return {"B": 0.0, "S": 0.0, "O": 0.0, "T": 0.0}


def percentages(beach, seabed, ocean, zone, weight, denom, n_emitted, n_present) -> dict:
    """Percent of emitted number, or of emitted mass when `weight` is mass."""

    def block(selected: np.ndarray) -> dict[str, float]:
        w = weight[selected]
        return {
            "B": 100.0 * w[beach[selected]].sum() / denom,
            "S": 100.0 * w[seabed[selected]].sum() / denom,
            "O": 100.0 * w[ocean[selected]].sum() / denom,
            "T": 100.0 * w.sum() / denom,
        }

    present = np.ones(len(weight), dtype=bool)
    total = block(present)
    total["OoD"] = 100.0 * (n_emitted - n_present) / n_emitted
    return {
        "total": total,
        "zones": {z: block(zone == z) for z in (1, 2, 3)},
    }


def summarise_case(case_dir: Path, stem: str, xml_path: Path, ring: np.ndarray) -> dict:
    files = tracer_files(case_dir, stem)
    if len(files) < 2:
        raise FileNotFoundError(f"{stem}: expected the initial and final VTU, found {len(files)}")
    initial = read_vtu(files[0])
    final = read_vtu(files[-1])
    n_releases = emission_count(xml_path)
    n_initial = len(initial["id"])
    n_emitted = n_initial * n_releases
    n_present = len(final["id"])
    if n_present > n_emitted:
        raise ValueError(f"{stem}: final VTU has {n_present} tracers but {n_emitted} were emitted")

    released_mass = float(particle_mass(initial["radius"], initial["density"]).sum())
    initial_mass = released_mass * n_releases
    mass = particle_mass(final["radius"], final["density"])
    points = final["Points"]
    beach, seabed, ocean = classify(points[:, 2], final["bathymetry"])
    zone = assign_zones(points[:, 0], points[:, 1], ring)
    number = percentages(
        beach, seabed, ocean, zone, np.ones(n_present), float(n_emitted), n_emitted, n_present
    )
    by_mass = percentages(
        beach, seabed, ocean, zone, mass, initial_mass, n_emitted, n_present
    )
    # Mass of tracers that left the domain is unknown. These runs keep them,
    # so OoD mass is zero whenever the final file still holds every release.
    if n_present != n_emitted:
        by_mass["total"]["OoD"] = float("nan")
    return {"number": number, "mass": by_mass, "n_present": n_present, "n_emitted": n_emitted}


def fnum(value: float) -> str:
    if not np.isfinite(value):
        return "{--}"
    return f"{value:.2f}"


def number_row(label: str, stats: dict) -> str:
    total = stats["total"]
    cells = [fnum(total[key]) for key in ("B", "S", "O", "OoD")]
    for zone_id in (1, 2, 3):
        zone = stats["zones"][zone_id]
        cells.extend(fnum(zone[key]) for key in ("B", "S", "O", "T"))
    return r"\textbf{" + label + "} & " + " & ".join(cells) + r" \\"


def mass_rows(label: str, case: dict) -> list[str]:
    rows = []
    for index, (kind, key) in enumerate((("mass", "mass"), ("number", "number"))):
        total = case[key]["total"]
        cells = [fnum(total[name]) for name in ("B", "S", "O", "OoD")]
        for zone_id in (1, 2, 3):
            zone = case[key]["zones"][zone_id]
            cells.extend(fnum(zone[name]) for name in ("B", "S", "O", "T"))
        body = " & ".join(cells) + r" \\"
        if index == 0:
            rows.append(
                r"\multirow{2}*{\textbf{" + label + r"}} & \textit{" + kind + "} & " + body
            )
        else:
            rows.append(r" & \textit{" + kind + "} & " + body)
    return rows


def table_header(column_spec: str, with_metric: bool) -> list[str]:
    if with_metric:
        zone_rules = r"\cmidrule(lr){3-6} \cmidrule(lr){7-10} \cmidrule(lr){11-14} \cmidrule(lr){15-18}"
        sub_rules = r"\cmidrule(lr){3-5} \cmidrule(lr){7-9} \cmidrule(lr){11-13} \cmidrule(lr){15-17}"
        titles = (
            r"& & \multicolumn{4}{c|}{\textbf{Total}} & "
            r"\multicolumn{4}{c|}{\textbf{Zona 1}} & "
            r"\multicolumn{4}{c|}{\textbf{Zona 2}} & "
            r"\multicolumn{4}{c}{\textbf{Zona 3}} \\"
        )
        names = (
            r"& & {B} & {S} & {O} & {OoD} & {B} & {S} & {O} & {\textbf{Total Z1}} & "
            r"{B} & {S} & {O} & {\textbf{Total Z2}} & {B} & {S} & {O} & {\textbf{Total Z3}} \\"
        )
    else:
        zone_rules = r"\cmidrule(lr){2-5} \cmidrule(lr){6-9} \cmidrule(lr){10-13} \cmidrule(lr){14-17}"
        sub_rules = r"\cmidrule(lr){6-8} \cmidrule(lr){10-12} \cmidrule(lr){14-16}"
        titles = (
            r"& \multicolumn{4}{c|}{\textbf{Total}} & "
            r"\multicolumn{4}{c|}{\textbf{Zona 1}} & "
            r"\multicolumn{4}{c|}{\textbf{Zona 2}} & "
            r"\multicolumn{4}{c}{\textbf{Zona 3}} \\"
        )
        names = (
            r"& {B} & {S} & {O} & {OoD} & {B} & {S} & {O} & {\textbf{Total Z1}} & "
            r"{B} & {S} & {O} & {\textbf{Total Z2}} & {B} & {S} & {O} & {\textbf{Total Z3}} \\"
        )
    return [
        r"\resizebox{\textwidth}{!}{%",
        r"\begin{tabular}{" + column_spec + "}",
        r"\toprule",
        titles,
        zone_rules,
        names,
        sub_rules,
        r"\midrule",
    ]


def table_footer() -> list[str]:
    return [
        r"\bottomrule",
        r"\end{tabular}%",
        r"}",
        r"\begin{itemize}%[noitemsep, topsep=5pt,itemsep=-0.0cm]",
        r"\tiny",
        r"    \item B: Beachbed, S: Seabed, O: Ocean, OoD: Out of Domain. ",
        r"    \item Z1: Zone 1, Z2: Zone 2, Z3: Zone 3.",
        r"\end{itemize}",
        r"\end{table}",
        "",
    ]


def open_table(
    caption: str, label: str, colored: bool = True, placement: str = "htbp"
) -> list[str]:
    lines = [
        r"\begin{table}[" + placement + "]",
        r"\centering",
    ]
    if colored:
        lines += [
            r"\ifdefined\AuthorTableColor",
            r"\AuthorTableColor",
            r"\fi",
        ]
    lines += [
        r"\caption{" + caption + "}",
        r"\label{" + label + "}",
    ]
    return lines


def padded_label(case_id: int) -> str:
    return f"S{case_id:02d}"


def fate_name(stats: dict) -> str:
    total = stats["total"]
    return "Beach" if total["B"] > total["S"] else "Seabed"


def zone_totals(stats: dict) -> list[str]:
    zones = stats["zones"]
    return [fnum(zones[zone_id]["T"]) for zone_id in (1, 2, 3)]


def compact_row(label: str, stats: dict, with_fate: bool = False) -> str:
    cells = zone_totals(stats)
    if with_fate:
        cells = [fate_name(stats), *cells]
    return r"\textbf{" + label + "} & " + " & ".join(cells) + r" \\"


def compact_mass_rows(label: str, case: dict) -> list[str]:
    rows = []
    for index, (kind, key) in enumerate((("mass", "mass"), ("number", "number"))):
        body = " & ".join(zone_totals(case[key])) + r" \\"
        if index == 0:
            rows.append(
                r"\multirow{2}*{\textbf{" + label + r"}} & \textit{" + kind + "} & " + body
            )
        else:
            rows.append(r" & \textit{" + kind + "} & " + body)
    return rows


def compact_table_header(with_fate: bool = False, with_metric: bool = False) -> list[str]:
    if with_metric:
        return [
            r"\begin{tabular}{l l S[table-format=3.2] S[table-format=3.2] S[table-format=3.2]}",
            r"\toprule",
            r"\textbf{Case} & & {Z1 \%} & {Z2 \%} & {Z3 \%} \\",
            r"\midrule",
        ]
    if with_fate:
        return [
            r"\begin{tabular}{l l S[table-format=3.2] S[table-format=3.2] S[table-format=3.2]}",
            r"\toprule",
            r"\textbf{Case} & \textbf{Fate} & {Z1 \%} & {Z2 \%} & {Z3 \%} \\",
            r"\midrule",
        ]
    return [
        r"\begin{tabular}{l S[table-format=3.2] S[table-format=3.2] S[table-format=3.2]}",
        r"\toprule",
        r"\textbf{Case} & {Z1 \%} & {Z2 \%} & {Z3 \%} \\",
        r"\midrule",
    ]


def compact_table_footer() -> list[str]:
    return [
        r"\bottomrule",
        r"\end{tabular}",
        r"",
        r"\vspace{4pt}",
        r"\noindent Z1: Zone 1, Z2: Zone 2, Z3: Zone 3.",
        r"\end{table}",
        "",
    ]


def grouped_compact_table(
    cases: list[dict],
    group: int,
    extra_rule_after: set[int] | None = None,
    with_fate: bool = False,
) -> list[str]:
    extra_rule_after = extra_rule_after or set()
    lines = []
    for index, case in enumerate(cases, start=1):
        lines.append(
            compact_row(
                padded_label(case["id"]), case["stats"]["number"], with_fate=with_fate
            )
        )
        if index != len(cases) and (index % group == 0 or index in extra_rule_after):
            lines.append(r"\midrule")
            if index in extra_rule_after:
                lines.append(r"\midrule")
    return lines


def document_preamble() -> list[str]:
    return [
        r"\documentclass{article}",
        r"\usepackage[margin=1.5cm]{geometry}",
        r"\usepackage{booktabs}",
        r"\usepackage{siunitx}",
        r"\usepackage{multirow}",
        r"\usepackage{graphicx}",
        r"\usepackage{xcolor}",
        r"\newif\ifshowrevisions",
        r"\showrevisionstrue",
        r"\newcommand{\AuthorTableColor}{\ifshowrevisions\color{orange}\fi}",
        r"\begin{document}",
        "",
    ]


def grouped_number_table(cases: list[tuple[str, dict]], group: int, extra_rule_after: set[int]) -> list[str]:
    lines = []
    for index, (label, case) in enumerate(cases, start=1):
        lines.append(number_row(label, case["number"]))
        if index != len(cases) and (index % group == 0 or index in extra_rule_after):
            lines.append(r"\midrule")
            if index in extra_rule_after:
                lines.append(r"\midrule")
    return lines


def discover_cases(root: Path) -> list[dict]:
    found = []
    for folder in sorted(root.glob("S??-*_out")):
        stem = folder.name[: -len("_out")]
        match = re.fullmatch(r"S(\d+)-(.+)", stem)
        if not match:
            continue
        found.append(
            {
                "id": int(match.group(1)),
                "kind": match.group(2),
                "label": f"S{int(match.group(1))}",
                "stem": stem,
                "folder": folder,
                "xml": root / f"{stem}.xml",
            }
        )
    return found


def build_document(cases: dict[int, dict]) -> str:
    material = [cases[i] for i in range(1, 13)]
    reference = [cases[i] for i in (5, 8, 11)]
    degradation = [cases[REFERENCE_ID], cases[13], cases[14]]
    method1 = [cases[i] for i in range(15, 33)]
    method2 = [cases[i] for i in range(33, 42)]
    diffusion = [cases[i] for i in (42, 43, 44, 45)]

    parts = document_preamble()

    parts += open_table(
        "Sensitivity analysis of MOHID-Lagrangian simulations illustrating the "
        "influence of tracer density and particle radius on transport behaviour, "
        "vertical redistribution, and deposition patterns under identical "
        "environmental forcing conditions.",
        "tab:Sensitivity analysis: density and radius effects.",
    )
    parts += table_header(NUMBER_SPEC, with_metric=False)
    parts += grouped_number_table(
        [(case["label"], case["stats"]) for case in material], group=3, extra_rule_after=set()
    )
    parts += table_footer()

    parts += open_table(
        "Sensitivity analysis of key particle properties in MOHID-Lagrangian, "
        "illustrating how variations in tracer degradation affect transport "
        "behavior, vertical redistribution, and deposition patterns under "
        "identical environmental forcing conditions.",
        "Sensitivity analysis: Degradation.",
    )
    parts += table_header(MASS_SPEC, with_metric=True)
    for index, case in enumerate(degradation):
        parts += mass_rows(case["label"], case["stats"])
        if index != len(degradation) - 1:
            parts.append(r"\midrule")
    parts += table_footer()

    parts += open_table(
        " Reference sensitivity-analysis results for scenarios S5, S8, and S11, "
        "extracted from the complete density- and radius-sensitivity analysis "
        "(Table \\ref{tab:Sensitivity analysis: density and radius effects.}) to "
        "provide a baseline for comparison with the resuspension simulations "
        "presented in the following Tables.",
        "tab:Sensitivity analysis: density and radius effects (S5,S8,S11).",
    )
    parts += table_header(NUMBER_SPEC, with_metric=False)
    for case in reference:
        parts.append(number_row(case["label"], case["stats"]["number"]))
    parts += table_footer()

    parts += open_table(
        "Sensitivity analysis of the resuspension parameterization using Method 1, "
        "showing the influence of sediment resuspension on tracer transport, "
        "redistribution, and deposition.",
        "Sensitivity analysis: Resuspension effect (Method 1).",
    )
    parts += table_header(NUMBER_SPEC, with_metric=False)
    parts += grouped_number_table(
        [(case["label"], case["stats"]) for case in method1],
        group=3,
        extra_rule_after={9},
    )
    parts += table_footer()

    parts += open_table(
        "Sensitivity analysis of the resuspension parameterization using Method 2, "
        "showing the influence of sediment resuspension on tracer transport, "
        "redistribution, and deposition.",
        "Sensitivity analysis: Resuspension effect (Method 2).",
    )
    parts += table_header(NUMBER_SPEC, with_metric=False)
    parts += grouped_number_table(
        [(case["label"], case["stats"]) for case in method2], group=3, extra_rule_after=set()
    )
    parts += table_footer()

    parts += open_table(
        "Sensitivity analysis of MOHID-Lagrangian simulations illustrating the "
        "influence of the horizontal and vertical diffusion coefficients on "
        "transport behaviour, vertical redistribution, and deposition patterns "
        "under identical environmental forcing conditions.",
        "Sensitivity analysis: Diffusion.",
    )
    parts += table_header(NUMBER_SPEC, with_metric=False)
    parts.append(number_row(cases[REFERENCE_ID]["label"], cases[REFERENCE_ID]["stats"]["number"]))
    parts.append(r"\midrule")
    parts.append(number_row(diffusion[0]["label"], diffusion[0]["stats"]["number"]))
    parts.append(number_row(diffusion[1]["label"], diffusion[1]["stats"]["number"]))
    parts.append(r"\midrule")
    parts.append(number_row(diffusion[2]["label"], diffusion[2]["stats"]["number"]))
    parts.append(number_row(diffusion[3]["label"], diffusion[3]["stats"]["number"]))
    parts += table_footer()

    parts.append(r"\end{document}")
    parts.append("")
    return "\n".join(parts)


def compact_preamble() -> list[str]:
    return [
        r"\documentclass{article}",
        r"\usepackage[margin=1.5cm]{geometry}",
        r"\usepackage{booktabs}",
        r"\usepackage{siunitx}",
        r"\usepackage{multirow}",
        r"\usepackage{float}",
        r"\begin{document}",
        "",
    ]


def compact_open(caption: str, label: str) -> list[str]:
    return open_table(caption, label, colored=False, placement="H")


def build_compact_document(cases: dict[int, dict]) -> str:
    material = [cases[i] for i in range(1, 13)]
    reference = [cases[i] for i in (5, 8, 11)]
    degradation = [cases[REFERENCE_ID], cases[13], cases[14]]
    method1 = [cases[i] for i in range(15, 33)]
    method2 = [cases[i] for i in range(33, 42)]
    diffusion = [cases[i] for i in (42, 43, 44, 45)]

    parts = compact_preamble()

    parts += [
        r"\subsubsection{Sensitivity analysis: density and radius effects}",
        "",
    ]
    parts += compact_open(
        "Sensitivity analysis of MOHID-Lagrangian simulations illustrating the "
        "influence of tracer density and particle radius on transport behaviour, "
        "vertical redistribution, and deposition patterns under identical "
        "environmental forcing conditions.",
        "tab:Sensitivity analysis: density and radius effects.",
    )
    parts += compact_table_header(with_fate=True)
    parts += grouped_compact_table(material, group=3, with_fate=True)
    parts += compact_table_footer()

    parts += [
        r"\subsubsection{Sensitivity analysis: degradation effects}",
        "",
    ]
    parts += compact_open(
        "Sensitivity analysis of key particle properties in MOHID-Lagrangian, "
        "illustrating how variations in tracer degradation affect transport "
        "behavior, vertical redistribution, and deposition patterns under "
        "identical environmental forcing conditions.",
        "Sensitivity analysis: Degradation.",
    )
    parts += compact_table_header(with_metric=True)
    for index, case in enumerate(degradation):
        parts += compact_mass_rows(padded_label(case["id"]), case["stats"])
        if index != len(degradation) - 1:
            parts.append(r"\midrule")
    parts += compact_table_footer()

    parts += [
        r"\subsubsection{Sensitivity analysis: resuspension effects}",
        "",
    ]
    parts += compact_open(
        " Reference sensitivity-analysis results for scenarios S5, S8, and S11, "
        "extracted from the complete density- and radius-sensitivity analysis "
        "(Table \\ref{tab:Sensitivity analysis: density and radius effects.}) to "
        "provide a baseline for comparison with the resuspension simulations "
        "presented in the following Tables.",
        "tab:Sensitivity analysis: density and radius effects (S5,S8,S11).",
    )
    parts += compact_table_header()
    for case in reference:
        parts.append(compact_row(padded_label(case["id"]), case["stats"]["number"]))
    parts += compact_table_footer()

    parts += [
        r"\paragraph{Method 1: Prescribed critical shear stress (S15--S32)}",
        "",
    ]
    parts += compact_open(
        "Sensitivity analysis of the resuspension parameterization using Method 1, "
        "showing the influence of sediment resuspension on tracer transport, "
        "redistribution, and deposition.",
        "Sensitivity analysis: Resuspension effect (Method 1).",
    )
    parts += compact_table_header()
    parts += grouped_compact_table(method1, group=3, extra_rule_after={9})
    parts += compact_table_footer()

    parts += [
        r"\paragraph{Method 2: Tracer-dependent critical shear stress (S33--S41)}",
        "",
    ]
    parts += compact_open(
        "Sensitivity analysis of the resuspension parameterization using Method 2, "
        "showing the influence of sediment resuspension on tracer transport, "
        "redistribution, and deposition.",
        "Sensitivity analysis: Resuspension effect (Method 2).",
    )
    parts += compact_table_header()
    parts += grouped_compact_table(method2, group=3)
    parts += compact_table_footer()

    parts += [
        r"\subsubsection{Sensitivity analysis: diffusion effects}",
        "",
    ]
    parts += compact_open(
        "Sensitivity analysis of MOHID-Lagrangian simulations illustrating the "
        "influence of the horizontal and vertical diffusion coefficients on "
        "transport behaviour, vertical redistribution, and deposition patterns "
        "under identical environmental forcing conditions.",
        "Sensitivity analysis: Diffusion.",
    )
    parts += compact_table_header()
    parts.append(compact_row(padded_label(REFERENCE_ID), cases[REFERENCE_ID]["stats"]["number"]))
    parts.append(r"\midrule")
    parts.append(compact_row(padded_label(diffusion[0]["id"]), diffusion[0]["stats"]["number"]))
    parts.append(compact_row(padded_label(diffusion[1]["id"]), diffusion[1]["stats"]["number"]))
    parts.append(r"\midrule")
    parts.append(compact_row(padded_label(diffusion[2]["id"]), diffusion[2]["stats"]["number"]))
    parts.append(compact_row(padded_label(diffusion[3]["id"]), diffusion[3]["stats"]["number"]))
    parts += compact_table_footer()

    parts.append(r"\end{document}")
    parts.append("")
    return "\n".join(parts)


def main() -> None:
    ring = load_polygon(CASE_DIR / "data" / "CambadosA1.kmz")
    discovered = discover_cases(CASE_DIR)
    if not discovered:
        raise SystemExit(f"No case outputs found in {CASE_DIR}")

    loaded: dict[int, dict] = {}
    for case in discovered:
        print(f"{case['label']} {case['kind']}")
        case["stats"] = summarise_case(case["folder"], case["stem"], case["xml"], ring)
        total = case["stats"]["number"]["total"]
        print(
            f"    n={case['stats']['n_present']}/{case['stats']['n_emitted']}  "
            f"B={total['B']:.2f} S={total['S']:.2f} O={total['O']:.2f} OoD={total['OoD']:.2f}"
        )
        loaded[case["id"]] = case

    expected = set(range(1, 46))
    missing = sorted(expected - set(loaded))
    if missing:
        names = ", ".join(f"S{i:02d}" for i in missing)
        raise SystemExit(f"Missing outputs for: {names}")

    OUT_TEX.write_text(build_document(loaded), encoding="utf-8")
    print(f"Wrote {OUT_TEX}")
    OUT_TEX_COMPACT.write_text(build_compact_document(loaded), encoding="utf-8")
    print(f"Wrote {OUT_TEX_COMPACT}")


if __name__ == "__main__":
    main()
