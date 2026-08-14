#!/usr/bin/env python3
"""Plot PLUME dispersion, polarization, resonance, and heating diagnostics."""

from __future__ import annotations

import argparse
import math
import re
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


REQUIRED_RES_KEYS = ("plume_input_file", "plume_output_dir", "figure_output_file")
BASE_COLUMNS = ("kperp", "kpar", "beta_ref", "vtpar_ref_over_c", "omega_r", "gamma")
CHANNEL_COLORS = {
    "total": "black",
    "LD": "#0072B2",
    "TTD": "#D55E00",
    "n=+1": "#009E73",
    "n=-1": "#CC79A7",
}


@dataclass
class ResConfig:
    plume_input_file: Path
    plume_output_dir: Path
    figure_output_file: Path
    selected_kperp: float | None = None
    selected_kpar: float | None = None
    selected_branch: int | None = None
    x_axis: str = "kpar"
    x_scale: str = "linear"
    normalize_velocity_to: str = "vthpar_species1"
    vpar_min: float = -4.0
    vpar_max: float = 4.0
    vperp_min: float = 0.0
    vperp_max: float = 4.0
    nvpar: int = 160
    nvperp: int = 100
    contour_quantity: str = "resonance_weight"
    use_log_abs_contours: bool = False
    vdf_contour_scale: str = "f"
    vdf_contours_per_decade: int = 3
    include_species: list[int] | str = "all"
    gamma_j_ylim: tuple[float, float] = (1.0e-9, 1.0)
    output_format: str | None = None
    dpi: int = 300
    show: bool = False
    raw: dict[str, Any] = field(default_factory=dict)


@dataclass
class SpeciesParams:
    index: int
    tau: float
    mu: float
    alpha: float
    charge_ratio: float
    density: float
    drift_va_ref: float

    @property
    def omega_over_ref(self) -> float:
        # PLUME uses muS=m_ref/m_s and Qs=q_ref/q_s, so Omega_s/Omega_ref=(q_s/m_s)/(q_ref/m_ref)=muS/Qs.
        return self.mu / self.charge_ratio

    @property
    def vthpar_over_ref(self) -> float:
        # vth_s/vth_ref = sqrt((T_s/T_ref) / (m_s/m_ref)) = sqrt(muS/tauS).
        return math.sqrt(self.mu / self.tau)


@dataclass
class PlumeInput:
    path: Path
    params: dict[str, Any]
    scans: list[dict[str, Any]]
    species: list[SpeciesParams]
    output_name: str | None
    data_name: str | None
    eigen: bool
    heating: bool
    low_n: bool


@dataclass
class OutputFiles:
    mode_files: dict[int, Path]
    dispersion_roots: Path | None = None
    dispersion_map: Path | None = None
    jet_plume_files: dict[tuple[int, int, str], Path] = field(default_factory=dict)


@dataclass
class DispersionData:
    path: Path
    data: np.ndarray
    columns: dict[str, int]
    branch: int
    eigen_available: bool = False
    heating_available: bool = False
    low_n_available: bool = False
    heat_start: int | None = None
    params_start: int | None = None


@dataclass
class Selection:
    index: int
    kperp: float
    kpar: float
    beta_ref: float
    omega_r: float
    gamma: float
    distance: float


@dataclass
class PolarizationDiagnostics:
    pxy: float = math.nan
    epar_fraction: float = math.nan


def parse_value(text: str) -> Any:
    value = text.strip().strip(",")
    low = value.lower()
    if low in {"true", "t", ".true.", "yes", "y"}:
        return True
    if low in {"false", "f", ".false.", "no", "n"}:
        return False
    if "," in value:
        return [parse_value(part) for part in value.split(",") if part.strip()]
    if (value.startswith("'") and value.endswith("'")) or (value.startswith('"') and value.endswith('"')):
        return value[1:-1]
    try:
        if any(c in value.lower() for c in (".", "e", "d")):
            return float(value.replace("D", "E").replace("d", "e"))
        return int(value)
    except ValueError:
        return value


def read_res_file(path: str | Path) -> ResConfig:
    res_path = Path(path).expanduser()
    values: dict[str, Any] = {}
    with res_path.open("r", encoding="utf-8") as handle:
        for lineno, line in enumerate(handle, start=1):
            clean = line.split("#", 1)[0].strip()
            if not clean:
                continue
            if "=" not in clean:
                raise ValueError(f"{res_path}:{lineno}: expected 'key = value', got {line.rstrip()!r}")
            key, value = clean.split("=", 1)
            values[key.strip().lower()] = parse_value(value)

    missing = [key for key in REQUIRED_RES_KEYS if key not in values]
    if missing:
        raise ValueError(f"{res_path}: missing required key(s): {', '.join(missing)}")

    base = res_path.parent

    def cfg_path(key: str) -> Path:
        p = Path(str(values[key])).expanduser()
        return p if p.is_absolute() else (base / p).resolve()

    include = values.get("include_species", "all")
    if isinstance(include, str) and include.lower() != "all":
        include = [int(x) for x in re.split(r"[,\s]+", include) if x]
    elif isinstance(include, list):
        include = [int(x) for x in include]

    x_axis = str(values.get("x_axis", "kpar")).lower()
    if x_axis not in {"kperp", "kpar", "k"}:
        raise ValueError("x_axis must be one of: kperp, kpar, k")
    x_scale = str(values.get("x_scale", "linear")).lower()
    if x_scale not in {"linear", "log"}:
        raise ValueError("x_scale must be one of: linear, log")
    gamma_j_ylim = _float_pair(values.get("gamma_j_ylim", [1.0e-9, 1.0]), "gamma_j_ylim")
    if gamma_j_ylim[0] <= 0.0 or gamma_j_ylim[1] <= gamma_j_ylim[0]:
        raise ValueError("gamma_j_ylim must be two positive increasing values, e.g. gamma_j_ylim = 1e-9, 1")
    vdf_contour_scale = str(values.get("vdf_contour_scale", "f")).lower()
    if vdf_contour_scale not in {"f", "log10_f"}:
        raise ValueError("vdf_contour_scale must be one of: f, log10_f")
    vdf_contours_per_decade = int(values.get("vdf_contours_per_decade", 3))
    if vdf_contours_per_decade < 1:
        raise ValueError("vdf_contours_per_decade must be >= 1")

    return ResConfig(
        plume_input_file=cfg_path("plume_input_file"),
        plume_output_dir=cfg_path("plume_output_dir"),
        figure_output_file=cfg_path("figure_output_file"),
        selected_kperp=_optional_float(values.get("selected_kperp")),
        selected_kpar=_optional_float(values.get("selected_kpar")),
        selected_branch=_optional_int(values.get("selected_branch")),
        x_axis=x_axis,
        x_scale=x_scale,
        normalize_velocity_to=str(values.get("normalize_velocity_to", "vthpar_species1")),
        vpar_min=float(values.get("vpar_min", -4.0)),
        vpar_max=float(values.get("vpar_max", 4.0)),
        vperp_min=float(values.get("vperp_min", 0.0)),
        vperp_max=float(values.get("vperp_max", 4.0)),
        nvpar=int(values.get("nvpar", 160)),
        nvperp=int(values.get("nvperp", 100)),
        contour_quantity=str(values.get("contour_quantity", "resonance_weight")),
        use_log_abs_contours=bool(values.get("use_log_abs_contours", False)),
        vdf_contour_scale=vdf_contour_scale,
        vdf_contours_per_decade=vdf_contours_per_decade,
        include_species=include,
        gamma_j_ylim=gamma_j_ylim,
        output_format=_optional_str(values.get("output_format")),
        dpi=int(values.get("dpi", 300)),
        show=bool(values.get("show", False)),
        raw=values,
    )


def _optional_float(value: Any) -> float | None:
    return None if value is None or str(value).strip() == "" else float(value)


def _optional_int(value: Any) -> int | None:
    return None if value is None or str(value).strip() == "" else int(value)


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _float_pair(value: Any, key: str) -> tuple[float, float]:
    if isinstance(value, list) and len(value) == 2:
        return (float(value[0]), float(value[1]))
    if isinstance(value, str):
        parts = [part for part in re.split(r"[,\s]+", value.strip()) if part]
        if len(parts) == 2:
            return (float(parts[0]), float(parts[1]))
    raise ValueError(f"{key} must contain two values, e.g. {key} = 1e-9, 1")


def read_plume_input(path: str | Path) -> PlumeInput:
    p = Path(path)
    namelists: dict[str, list[dict[str, Any]]] = {}
    current_name: str | None = None
    current_values: dict[str, Any] = {}

    for raw_line in p.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.split("!", 1)[0].strip()
        if not line:
            continue
        if line.startswith("&"):
            current_name = line[1:].strip().lower()
            current_values = {}
            continue
        if line.startswith("/"):
            if current_name is not None:
                namelists.setdefault(current_name, []).append(current_values)
            current_name = None
            current_values = {}
            continue
        if current_name is None:
            continue
        for chunk in _split_assignments(line):
            if "=" not in chunk:
                continue
            key, value = chunk.split("=", 1)
            current_values[key.strip().lower()] = parse_value(value)

    params = namelists.get("params", [{}])[0]
    nspec = int(params.get("nspec", 0))
    if nspec < 1:
        raise ValueError(f"{p}: could not determine nspec from &params")

    species: list[SpeciesParams] = []
    for idx in range(1, nspec + 1):
        block = namelists.get(f"species_{idx}", [{}])[0]
        try:
            species.append(
                SpeciesParams(
                    index=idx,
                    tau=float(block["taus"]),
                    mu=float(block["mus"]),
                    alpha=float(block["alphs"]),
                    charge_ratio=float(block["qs"]),
                    density=float(block["ds"]),
                    drift_va_ref=float(block.get("vvs", 0.0)),
                )
            )
        except KeyError as exc:
            raise ValueError(f"{p}: species_{idx} is missing required key {exc.args[0]!r}") from exc

    scans = []
    for key, vals in sorted(namelists.items()):
        if key.startswith("scan_input_"):
            scans.extend(vals)

    eigen = any(bool(scan.get("eigen", False)) for scan in scans)
    heating = any(bool(scan.get("heating", scan.get("heat", False))) for scan in scans)
    low_n = bool(params.get("low_n", True))
    return PlumeInput(
        path=p,
        params=params,
        scans=scans,
        species=species,
        output_name=_optional_str(params.get("outputname")),
        data_name=_optional_str(params.get("dataname")),
        eigen=eigen,
        heating=heating,
        low_n=low_n,
    )


def _split_assignments(line: str) -> list[str]:
    parts: list[str] = []
    current: list[str] = []
    quote: str | None = None
    for ch in line:
        if ch in {"'", '"'}:
            quote = None if quote == ch else ch if quote is None else quote
        if ch == "," and quote is None:
            parts.append("".join(current).strip())
            current = []
        else:
            current.append(ch)
    if current:
        parts.append("".join(current).strip())
    return parts


def discover_plume_outputs(config: ResConfig, plume_input: PlumeInput) -> OutputFiles:
    output_path = config.plume_output_dir
    outdir = output_path.parent if output_path.is_file() else output_path
    stem = plume_input.output_name or config.plume_input_file.stem
    mode_files: dict[int, Path] = {}

    if output_path.is_file():
        match = re.search(r"\.mode(\d+)$", output_path.name)
        if not match:
            raise ValueError(f"plume_output_dir points to a file, but it is not a .modeN file: {output_path}")
        mode_files[int(match.group(1))] = output_path
    else:
        for path in sorted(outdir.glob("*.mode*")):
            if path.name.startswith("."):
                continue
            match = re.search(r"\.mode(\d+)$", path.name)
            if not match:
                continue
            if stem and stem not in path.name:
                continue
            mode_files.setdefault(int(match.group(1)), path)
        if not mode_files:
            for path in sorted(outdir.glob("*.mode*")):
                if path.name.startswith("."):
                    continue
                match = re.search(r"\.mode(\d+)$", path.name)
                if match:
                    mode_files.setdefault(int(match.group(1)), path)

    roots = _first_existing(sorted(outdir.glob(f"dispersion_{stem}.roots"))) or _first_existing(sorted(outdir.glob("dispersion_*.roots")))
    dmap = _first_existing(sorted(outdir.glob(f"dispersion_{stem}.map"))) or _first_existing(sorted(outdir.glob("dispersion_*.map")))

    jet_files: dict[tuple[int, int, str], Path] = {}
    for path in outdir.glob("*"):
        m = re.search(r"\.(cpar|cperp|cparcart|cperp1|cperp2|df1gyro\.real|df1gyro\.imag|df1\.real|df1\.imag)\.specie(\d+)\.mode(\d+)$", path.name)
        if m:
            quantity, species, mode = m.groups()
            jet_files[(int(species), int(mode), quantity)] = path

    return OutputFiles(mode_files=mode_files, dispersion_roots=roots, dispersion_map=dmap, jet_plume_files=jet_files)


def _first_existing(paths: list[Path]) -> Path | None:
    return paths[0] if paths else None


def read_dispersion_output(path: Path, plume_input: PlumeInput, branch: int) -> DispersionData:
    data = np.loadtxt(path, ndmin=2)
    if data.shape[1] < 6:
        raise ValueError(f"{path}: expected at least 6 columns, found {data.shape[1]}")

    columns = {name: i for i, name in enumerate(BASE_COLUMNS)}
    nspec = len(plume_input.species)
    eigen_block_cols = 12 + 8 * nspec
    eigen_available = plume_input.eigen and data.shape[1] >= 6 + 12
    if eigen_available:
        for name, col in {
            "Bx_re": 6,
            "Bx_im": 7,
            "By_re": 8,
            "By_im": 9,
            "Bz_re": 10,
            "Bz_im": 11,
            "Ex_re": 12,
            "Ex_im": 13,
            "Ey_re": 14,
            "Ey_im": 15,
            "Ez_re": 16,
            "Ez_im": 17,
        }.items():
            if col < data.shape[1]:
                columns[name] = col

    # Per output.md, P_j starts immediately after the full eigen block:
    # 1-based P_j column = 18 + 8*nspec + j, so 0-based start = 18 + 8*nspec.
    heat_start = 6 + eigen_block_cols if eigen_available else 6
    heating_available = data.shape[1] >= heat_start + nspec
    split_start = heat_start + nspec
    low_n_available = data.shape[1] >= split_start + 7 * nspec
    params_start = split_start + 7 * nspec if low_n_available and data.shape[1] >= split_start + 7 * nspec + 6 * nspec else None

    if heating_available:
        for j in range(nspec):
            columns[f"P_species_{j+1}"] = heat_start + j
        if low_n_available:
            for j in range(nspec):
                offset = split_start + 7 * j
                labels = ("TTD_yy", "TTD_yz", "LD_zy", "LD_zz", "n0", "n+1", "n-1")
                for label_i, label in enumerate(labels):
                    columns[f"P_species_{j+1}_{label}"] = offset + label_i

    return DispersionData(
        path=path,
        data=data,
        columns=columns,
        branch=branch,
        eigen_available=eigen_available,
        heating_available=heating_available,
        low_n_available=low_n_available,
        heat_start=heat_start if heating_available else None,
        params_start=params_start,
    )


def read_species_resonance_outputs(files: OutputFiles, branch: int) -> dict[tuple[int, str], np.ndarray]:
    loaded: dict[tuple[int, str], np.ndarray] = {}
    for (species, mode, quantity), path in files.jet_plume_files.items():
        if mode != branch:
            continue
        try:
            loaded[(species, quantity)] = np.loadtxt(path, ndmin=2)
        except Exception as exc:  # pragma: no cover - diagnostic only
            warnings.warn(f"Could not read JET-PLUME file {path}: {exc}", RuntimeWarning)
    return loaded


def select_k_point(config: ResConfig, disp: DispersionData) -> Selection:
    data = disp.data
    kperp = data[:, disp.columns["kperp"]]
    kpar = data[:, disp.columns["kpar"]]
    target_kperp = config.selected_kperp if config.selected_kperp is not None else float(kperp[0])
    target_kpar = config.selected_kpar if config.selected_kpar is not None else float(kpar[0])
    kperp_scale = max(np.nanmax(kperp) - np.nanmin(kperp), np.nanmax(np.abs(kperp)), 1.0)
    kpar_scale = max(np.nanmax(kpar) - np.nanmin(kpar), np.nanmax(np.abs(kpar)), 1.0)
    distance = ((kperp - target_kperp) / kperp_scale) ** 2 + ((kpar - target_kpar) / kpar_scale) ** 2
    idx = int(np.nanargmin(distance))
    omega = float(data[idx, disp.columns["omega_r"]])
    gamma = float(data[idx, disp.columns["gamma"]])
    if not np.isfinite([omega, gamma]).all():
        raise ValueError(f"Selected row {idx} has non-finite omega/gamma")
    beta_ref = float(data[idx, disp.columns["beta_ref"]])
    return Selection(idx, float(kperp[idx]), float(kpar[idx]), beta_ref, omega, gamma, float(math.sqrt(distance[idx])))


def compute_species_velocity_grids(config: ResConfig) -> tuple[np.ndarray, np.ndarray]:
    vpar = np.linspace(config.vpar_min, config.vpar_max, config.nvpar)
    vperp = np.linspace(config.vperp_min, config.vperp_max, config.nvperp)
    return np.meshgrid(vpar, vperp)


def velocity_normalization(config: ResConfig, plume_input: PlumeInput) -> float:
    token = config.normalize_velocity_to.lower()
    if token in {"vthpar_species1", "vth_species1", "vth1"}:
        return plume_input.species[0].vthpar_over_ref
    match = re.fullmatch(r"vthpar_species(\d+)", token)
    if match:
        idx = int(match.group(1))
        if idx < 1 or idx > len(plume_input.species):
            raise ValueError(f"normalize_velocity_to requested species {idx}, but only {len(plume_input.species)} exist")
        return plume_input.species[idx - 1].vthpar_over_ref
    if token in {"vthpar_ref", "vth_ref", "ref"}:
        return 1.0
    raise ValueError(f"Unsupported normalize_velocity_to={config.normalize_velocity_to!r}")


def fallback_vdf_values(
    species: SpeciesParams,
    selection: Selection,
    norm: float,
    vpar_grid: np.ndarray,
    vperp_grid: np.ndarray,
) -> np.ndarray:
    vpar_width = max(species.vthpar_over_ref / norm, 1.0e-30)
    vperp_width = max(math.sqrt(species.alpha) * species.vthpar_over_ref / norm, 1.0e-30)
    drift = 0.0
    if selection.beta_ref > 0.0:
        # vvS is normalized to v_A,ref and beta_ref=vth_ref^2/v_A,ref^2.
        drift = species.drift_va_ref / (math.sqrt(selection.beta_ref) * norm)
    return species.density * np.exp(-((vpar_grid - drift) / vpar_width) ** 2 - (vperp_grid / vperp_width) ** 2)


def fallback_vdf_contour_levels(
    config: ResConfig,
    species_list: list[SpeciesParams],
    selection: Selection,
    norm: float,
    vpar_grid: np.ndarray,
    vperp_grid: np.ndarray,
) -> np.ndarray:
    values = [
        fallback_vdf_values(species, selection, norm, vpar_grid, vperp_grid)
        for species in species_list
    ]
    positive = np.concatenate([grid[np.isfinite(grid) & (grid > 0.0)] for grid in values])
    if positive.size == 0:
        return np.linspace(0.0, 1.0, 16)
    min_decade = math.floor(math.log10(float(np.min(positive))))
    max_decade = max(0, math.ceil(math.log10(float(np.max(positive)))))
    log_levels = np.arange(
        min_decade,
        max_decade + 1.0 / config.vdf_contours_per_decade,
        1.0 / config.vdf_contours_per_decade,
    )
    if config.vdf_contour_scale == "log10_f":
        return log_levels
    return 10.0**log_levels


def compute_resonance_curves_or_weights(
    species: SpeciesParams,
    selection: Selection,
    norm_vth_over_ref: float,
) -> dict[int, float]:
    if abs(selection.kpar) < 1.0e-12:
        warnings.warn("Selected k_parallel is too close to zero; resonance velocities are undefined.", RuntimeWarning)
        return {}
    curves: dict[int, float] = {}
    for n in (0, 1, -1):
        v_ref_units = (selection.omega_r - n * species.omega_over_ref) / selection.kpar
        curves[n] = v_ref_units / norm_vth_over_ref
    return curves


def x_values(disp: DispersionData, x_axis: str) -> np.ndarray:
    kperp = disp.data[:, disp.columns["kperp"]]
    kpar = disp.data[:, disp.columns["kpar"]]
    if x_axis == "kperp":
        return kperp
    if x_axis == "kpar":
        return kpar
    return np.sqrt(kperp**2 + kpar**2)


def plot_signed_log(ax: plt.Axes, x: np.ndarray, y: np.ndarray, color: str, label: str, linewidth: float = 1.7) -> None:
    y = np.asarray(y, dtype=float)
    finite = np.isfinite(x) & np.isfinite(y) & (y != 0.0)
    if not np.any(finite):
        return
    xf = x[finite]
    yf = y[finite]
    start = 0
    signs = np.sign(yf)
    for i in range(1, len(yf) + 1):
        if i == len(yf) or signs[i] != signs[start]:
            if i - start > 0:
                ax.plot(
                    xf[start:i],
                    np.abs(yf[start:i]),
                    color=color,
                    linestyle="-" if signs[start] > 0 else "--",
                    linewidth=linewidth,
                    label=label if start == 0 else None,
                )
            start = i


def selected_polarization_diagnostics(disp: DispersionData, selection: Selection) -> PolarizationDiagnostics:
    required = ("Ex_re", "Ex_im", "Ey_re", "Ey_im", "Ez_re", "Ez_im")
    if not (disp.eigen_available and all(name in disp.columns for name in required)):
        return PolarizationDiagnostics()
    row = disp.data[selection.index]
    ex = row[disp.columns["Ex_re"]] + 1j * row[disp.columns["Ex_im"]]
    ey = row[disp.columns["Ey_re"]] + 1j * row[disp.columns["Ey_im"]]
    ez = row[disp.columns["Ez_re"]] + 1j * row[disp.columns["Ez_im"]]
    denom_xy = abs(ex) ** 2 + abs(ey) ** 2
    denom_total = math.sqrt(abs(ex) ** 2 + abs(ey) ** 2 + abs(ez) ** 2)
    pxy = 2.0 * np.imag(ex * np.conj(ey)) / denom_xy if denom_xy > 0.0 else math.nan
    epar_fraction = abs(ez) / denom_total if denom_total > 0.0 else math.nan
    return PolarizationDiagnostics(float(pxy), float(epar_fraction))


def make_resonance_figure(
    config: ResConfig,
    plume_input: PlumeInput,
    disp: DispersionData,
    selection: Selection,
    jet_data: dict[tuple[int, str], np.ndarray],
) -> None:
    species_indices = resolve_species_indices(config, plume_input)
    nrows = 1 + len(species_indices)
    fig = plt.figure(figsize=(10.8, 3.0 + 2.4 * len(species_indices)), constrained_layout=True)
    gs = fig.add_gridspec(nrows=nrows, ncols=2, height_ratios=[1.0] + [1.05] * len(species_indices))

    x = x_values(disp, config.x_axis)
    selected_x = x[selection.index]
    pol_diag = selected_polarization_diagnostics(disp, selection)

    ax_disp = fig.add_subplot(gs[0, 0])
    plot_signed_log(ax_disp, x, disp.data[:, disp.columns["omega_r"]], "black", r"$\omega_r/\Omega_{\rm ref}$")
    plot_signed_log(ax_disp, x, disp.data[:, disp.columns["gamma"]], "red", r"$\gamma/\Omega_{\rm ref}$")
    ax_disp.axvline(selected_x, color="0.55", linewidth=0.8)
    ax_disp.set_yscale("log")
    apply_x_scale(ax_disp, x, config.x_scale)
    ax_disp.set_xlabel(x_label(config.x_axis))
    ax_disp.set_ylabel("absolute value")
    ax_disp.grid(True, which="both", linewidth=0.35, alpha=0.35)
    ax_disp.legend(fontsize=8, frameon=False)

    ax_pol = fig.add_subplot(gs[0, 1])
    if disp.eigen_available and all(name in disp.columns for name in ("Ex_re", "Ex_im", "Ey_re", "Ey_im", "Ez_re", "Ez_im")):
        ex = disp.data[:, disp.columns["Ex_re"]] + 1j * disp.data[:, disp.columns["Ex_im"]]
        ey = disp.data[:, disp.columns["Ey_re"]] + 1j * disp.data[:, disp.columns["Ey_im"]]
        ez = disp.data[:, disp.columns["Ez_re"]] + 1j * disp.data[:, disp.columns["Ez_im"]]
        denom = np.abs(ex) ** 2 + np.abs(ey) ** 2
        pxy = np.divide(2.0 * np.imag(ex * np.conj(ey)), denom, out=np.full_like(denom, np.nan), where=denom > 0)
        efrac = np.divide(np.abs(ez), np.sqrt(np.abs(ex) ** 2 + np.abs(ey) ** 2 + np.abs(ez) ** 2), out=np.full_like(denom, np.nan), where=(np.abs(ex) + np.abs(ey) + np.abs(ez)) > 0)
        pol_line = ax_pol.plot(x, pxy, color="#0072B2", linewidth=1.5, label=r"$\mathcal{P}_{xy}$")
        ax_pol.set_ylim(-1.05, 1.05)
        ax_pol.set_ylabel(r"$\mathcal{P}_{xy}$")
        ax2 = ax_pol.twinx()
        epar_line = ax2.plot(x, efrac, color="#D55E00", linewidth=1.4, label=r"$|E_\parallel|/|E|$")
        ax2.set_yscale("log")
        ax2.set_ylabel(r"$|E_\parallel|/|E|$")
        ax2.axvline(selected_x, color="0.55", linewidth=0.8)
        ax_pol.legend(handles=pol_line + epar_line, fontsize=8, frameon=False, loc="best")
    else:
        ax_pol.text(0.5, 0.5, "Eigenfunction columns unavailable", ha="center", va="center", transform=ax_pol.transAxes)
    ax_pol.axvline(selected_x, color="0.55", linewidth=0.8)
    apply_x_scale(ax_pol, x, config.x_scale)
    ax_pol.set_xlabel(x_label(config.x_axis))
    ax_pol.grid(True, linewidth=0.35, alpha=0.35)

    norm = velocity_normalization(config, plume_input)
    vpar_grid, vperp_grid = compute_species_velocity_grids(config)
    selected_species = [plume_input.species[idx - 1] for idx in species_indices]
    vdf_levels = fallback_vdf_contour_levels(config, selected_species, selection, norm, vpar_grid, vperp_grid)
    for row, species_idx in enumerate(species_indices, start=1):
        species = plume_input.species[species_idx - 1]
        ax_v = fig.add_subplot(gs[row, 0])
        draw_velocity_panel(ax_v, config, species, selection, norm, vpar_grid, vperp_grid, vdf_levels, pol_diag, jet_data)
        ax_g = fig.add_subplot(gs[row, 1], sharex=ax_disp)
        draw_gamma_panel(ax_g, config, disp, species, x, selected_x)

    fig.suptitle(
        f"{config.plume_input_file.stem} mode {disp.branch}: "
        rf"$k_\perp\rho_{{ref}}={selection.kperp:.3g}$, "
        rf"$k_\parallel\rho_{{ref}}={selection.kpar:.3g}$",
        fontsize=11,
    )
    config.figure_output_file.parent.mkdir(parents=True, exist_ok=True)
    save_kwargs = {"dpi": config.dpi}
    if config.output_format:
        save_kwargs["format"] = config.output_format
    fig.savefig(config.figure_output_file, **save_kwargs)
    if config.show:
        plt.show()
    plt.close(fig)


def resolve_species_indices(config: ResConfig, plume_input: PlumeInput) -> list[int]:
    if config.include_species == "all":
        return [s.index for s in plume_input.species]
    indices = list(config.include_species)
    bad = [idx for idx in indices if idx < 1 or idx > len(plume_input.species)]
    if bad:
        raise ValueError(f"Requested species {bad}, but available species are 1..{len(plume_input.species)}")
    return indices


def x_label(x_axis: str) -> str:
    labels = {
        "kperp": r"$k_\perp\rho_{\rm ref}$",
        "kpar": r"$k_\parallel\rho_{\rm ref}$",
        "k": r"$k\rho_{\rm ref}$",
    }
    return labels[x_axis]


def apply_x_scale(ax: plt.Axes, x: np.ndarray, scale: str) -> None:
    if scale == "log":
        if np.any(np.isfinite(x) & (x <= 0.0)):
            warnings.warn("x_scale=log requested, but some x values are <= 0 and will not be visible.", RuntimeWarning)
        ax.set_xscale("log")


def resonance_line_style(n: int, species: SpeciesParams, pol_diag: PolarizationDiagnostics) -> str:
    if n == 0:
        return "-" if np.isfinite(pol_diag.epar_fraction) and pol_diag.epar_fraction > 1.0e-4 else "--"

    pxy = pol_diag.pxy
    if not np.isfinite(pxy) or species.charge_ratio == 0.0:
        return "--"
    charge_sign = np.sign(species.charge_ratio)
    if n == 1:
        active = (pxy > 0.9 and charge_sign > 0.0) or (pxy < -0.9 and charge_sign < 0.0)
    elif n == -1:
        active = (pxy > 0.0 and charge_sign < 0.0) or (pxy < -0.9 and charge_sign > 0.0)
    else:
        active = False
    return "-" if active else "--"


def draw_velocity_panel(
    ax: plt.Axes,
    config: ResConfig,
    species: SpeciesParams,
    selection: Selection,
    norm: float,
    vpar_grid: np.ndarray,
    vperp_grid: np.ndarray,
    vdf_levels: np.ndarray,
    pol_diag: PolarizationDiagnostics,
    jet_data: dict[tuple[int, str], np.ndarray],
) -> None:
    quantity_map = {
        "gamma_density": "cpar",
        "field_particle_signature": "cpar",
        "resonance_weight": "cpar",
    }
    jet_key = quantity_map.get(config.contour_quantity, config.contour_quantity)
    grid = jet_data.get((species.index, jet_key))
    if grid is not None:
        plot_grid = np.abs(grid) if config.use_log_abs_contours else grid
        gvpar = np.linspace(config.vpar_min, config.vpar_max, plot_grid.shape[1])
        gvperp = np.linspace(config.vperp_min, config.vperp_max, plot_grid.shape[0])
        gvpar_grid, gvperp_grid = np.meshgrid(gvpar, gvperp)
        ax.contourf(gvpar_grid, gvperp_grid, plot_grid, levels=24, cmap="RdBu_r")
    else:
        # No JET-PLUME velocity-space output was discovered; show a light
        # density-weighted bi-Maxwellian reference in the requested velocity units.
        weight = fallback_vdf_values(species, selection, norm, vpar_grid, vperp_grid)
        plot_weight = np.log10(np.clip(weight, 1.0e-300, None)) if config.vdf_contour_scale == "log10_f" else weight
        ax.contourf(vpar_grid, vperp_grid, plot_weight, levels=vdf_levels, cmap="Greys", alpha=0.35)
        ax.text(0.02, 0.96, "analytic resonance overlay", transform=ax.transAxes, ha="left", va="top", fontsize=8)

    resonances = compute_resonance_curves_or_weights(species, selection, norm)
    resonance_styles = {
        0: ("LD/TTD n=0", CHANNEL_COLORS["LD"]),
        1: ("n=+1", CHANNEL_COLORS["n=+1"]),
        -1: ("n=-1", CHANNEL_COLORS["n=-1"]),
    }
    for n, vres in resonances.items():
        if config.vpar_min <= vres <= config.vpar_max:
            label, color = resonance_styles[n]
            linestyle = resonance_line_style(n, species, pol_diag)
            ax.axvline(vres, color=color, linestyle=linestyle, linewidth=1.5, label=label)
    ax.set_xlim(config.vpar_min, config.vpar_max)
    ax.set_ylim(config.vperp_min, config.vperp_max)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"$v_\parallel/v_{\rm norm}$")
    ax.set_ylabel(r"$v_\perp/v_{\rm norm}$")
    ax.set_title(
        rf"species {species.index}: $\Omega_s/\Omega_{{ref}}={species.omega_over_ref:.3g}$, "
        rf"$v_{{th,\parallel s}}/v_{{th,ref}}={species.vthpar_over_ref:.3g}$",
        fontsize=9,
    )
    ax.legend(fontsize=7, frameon=False, loc="upper right")


def draw_gamma_panel(ax: plt.Axes, config: ResConfig, disp: DispersionData, species: SpeciesParams, x: np.ndarray, selected_x: float) -> None:
    omega = disp.data[:, disp.columns["omega_r"]]
    gamma = disp.data[:, disp.columns["gamma"]]
    sign_omega = np.sign(np.where(omega == 0.0, np.nan, omega))
    with np.errstate(divide="ignore", invalid="ignore"):
        modal_gamma = (gamma / omega) * sign_omega
    plot_signed_log(ax, x, modal_gamma, "0.55", r"$\gamma/\omega_r$", linewidth=1.2)

    if disp.heating_available:
        p_total_col = disp.columns.get(f"P_species_{species.index}")
        if p_total_col is not None:
            values = -disp.data[:, p_total_col] / sign_omega
            plot_signed_log(ax, x, values, CHANNEL_COLORS["total"], rf"$\gamma_{{{species.index}}}/\omega_r$", linewidth=1.5)
        if disp.low_n_available:
            channel_cols = {
                "LD": [f"P_species_{species.index}_LD_zy", f"P_species_{species.index}_LD_zz"],
                "TTD": [f"P_species_{species.index}_TTD_yy", f"P_species_{species.index}_TTD_yz"],
                "n=+1": [f"P_species_{species.index}_n+1"],
                "n=-1": [f"P_species_{species.index}_n-1"],
            }
            channel_values: dict[str, np.ndarray] = {}
            for label, names in channel_cols.items():
                cols = [disp.columns[name] for name in names if name in disp.columns]
                if cols:
                    channel_values[label] = -np.nansum(disp.data[:, cols], axis=1) / sign_omega
            for label, values in channel_values.items():
                plot_signed_log(ax, x, values, CHANNEL_COLORS[label], label, linewidth=1.2)
        else:
            ax.text(0.02, 0.92, "low-n channels unavailable", transform=ax.transAxes, ha="left", va="top", fontsize=8)
    else:
        ax.text(0.02, 0.92, "species heating columns unavailable", transform=ax.transAxes, ha="left", va="top", fontsize=8)

    ax.axvline(selected_x, color="0.55", linewidth=0.8)
    ax.set_yscale("log")
    ax.set_ylim(*config.gamma_j_ylim)
    apply_x_scale(ax, x, config.x_scale)
    ax.set_xlabel(x_label(config.x_axis))
    ax.set_ylabel("absolute value; dashed is negative")
    ax.grid(True, which="both", linewidth=0.35, alpha=0.35)
    ax.legend(fontsize=7, frameon=False, ncol=2)


def dry_run_report(config: ResConfig, plume_input: PlumeInput, outputs: OutputFiles, disp: DispersionData, selection: Selection, jet_data: dict[tuple[int, str], np.ndarray]) -> None:
    print("Parsed .resplot config:")
    for key, value in config.raw.items():
        print(f"  {key} = {value}")
    print("\nPLUME input:")
    print(f"  file = {plume_input.path}")
    print(f"  outputName = {plume_input.output_name}")
    print(f"  dataName = {plume_input.data_name}")
    print(f"  nspec = {len(plume_input.species)}")
    print(f"  eigen requested = {plume_input.eigen}")
    print(f"  heating requested = {plume_input.heating}")
    print(f"  low_n = {plume_input.low_n}")
    for species in plume_input.species:
        print(
            f"  species {species.index}: tau={species.tau:g}, mu={species.mu:g}, alpha={species.alpha:g}, "
            f"Qs={species.charge_ratio:g}, Ds={species.density:g}, vvS={species.drift_va_ref:g}, "
            f"Omega_s/Omega_ref={species.omega_over_ref:g}, vth_s/vth_ref={species.vthpar_over_ref:g}"
        )

    print("\nDetected output files:")
    print(f"  roots = {outputs.dispersion_roots}")
    print(f"  map = {outputs.dispersion_map}")
    for branch, path in sorted(outputs.mode_files.items()):
        print(f"  mode {branch} = {path}")
    print(f"  JET-PLUME velocity files = {len(outputs.jet_plume_files)}")

    print("\nSelected dispersion file:")
    print(f"  branch = {disp.branch}")
    print(f"  file = {disp.path}")
    print(f"  shape = {disp.data.shape}")
    print(f"  available columns = {', '.join(sorted(disp.columns))}")
    print(f"  eigen block detected = {disp.eigen_available}")
    print(f"  heating block detected = {disp.heating_available}")
    print(f"  low-n channel block detected = {disp.low_n_available}")
    if disp.heating_available:
        p_columns = [
            f"species {idx}: {disp.columns[f'P_species_{idx}'] + 1}"
            for idx in range(1, len(plume_input.species) + 1)
            if f"P_species_{idx}" in disp.columns
        ]
        print(f"  P_j columns (1-based) = {', '.join(p_columns)}")

    print("\nSelected k point:")
    print(f"  row index = {selection.index}")
    print(f"  kperp = {selection.kperp:.8g}")
    print(f"  kpar = {selection.kpar:.8g}")
    print(f"  beta_ref = {selection.beta_ref:.8g}")
    print(f"  omega_r = {selection.omega_r:.8g}")
    print(f"  gamma = {selection.gamma:.8g}")
    print(f"  normalized nearest-neighbor distance = {selection.distance:.3e}")

    channels = ["total species P_j"] if disp.heating_available else []
    if disp.low_n_available:
        channels.extend(["LD", "TTD", "n=+1", "n=-1"])
    print(f"\nAvailable species/resonance channels: {', '.join(channels) if channels else 'none'}")
    if not jet_data:
        print("Velocity-space FPC contours: none found; figure will use analytic resonance overlays.")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("resplot", help="Path to a .resplot file")
    parser.add_argument("--dry-run", action="store_true", help="Parse and report available data without making the figure")
    args = parser.parse_args(argv)

    config = read_res_file(args.resplot)
    plume_input = read_plume_input(config.plume_input_file)
    outputs = discover_plume_outputs(config, plume_input)
    if not outputs.mode_files:
        raise SystemExit(f"No .modeN files found in {config.plume_output_dir}")
    branch = config.selected_branch or sorted(outputs.mode_files)[0]
    if branch not in outputs.mode_files:
        raise SystemExit(f"selected_branch={branch} is unavailable; found branches {sorted(outputs.mode_files)}")
    disp = read_dispersion_output(outputs.mode_files[branch], plume_input, branch)
    selection = select_k_point(config, disp)
    jet_data = read_species_resonance_outputs(outputs, branch)

    if args.dry_run:
        dry_run_report(config, plume_input, outputs, disp, selection, jet_data)
        return 0

    if not disp.low_n_available:
        warnings.warn("Low-n resonance decomposition columns were not detected; channel curves will be omitted.", RuntimeWarning)
    if not jet_data:
        warnings.warn("No JET-PLUME velocity-space contour files detected; using analytic resonance overlays only.", RuntimeWarning)
    make_resonance_figure(config, plume_input, disp, selection, jet_data)
    print(f"Saved {config.figure_output_file}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2)
