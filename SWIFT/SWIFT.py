#!/usr/bin/env python3
"""
plume_wave_viz.py

Single-file utilities to:
  - Parse a single eigenmode row from a PLUME ASCII .modeN file (no header).
  - Reconstruct the real-space magnetic fluctuation field \delta B(x,t) from the
    complex eigenvector and complex frequency.
  - Reconstruct the reference-species velocity fluctuation \delta U_ref(x,t).
  - Build a simple Cartesian sampling domain aligned with B0 = B0 z-hat.
  - Visualize static and animated 3D magnetic-field and velocity fluctuations with PyVista.

PyVista is optional: the core parsing/reconstruction/domain helpers work without it,
and visualization functions raise a clear error only if they are called.
This script is intended to be readable and easy to modify.

Primary format & convention references (as requested):
  - PLUME output column ordering (.modeN files):
      https://github.com/kgklein/PLUME/blob/main/output.md
    In particular, for .modeN from om_scan:
      cols 1-6: k_perp*rho_ref, k_par*rho_ref, beta_ref_par, vt_ref_par/c, omega_r/Omega_ref, gamma/Omega_ref
      if eigen=true:
        cols 7-12:   Re[Bx], Im[Bx], Re[By], Im[By], Re[Bz], Im[Bz]
        cols 13-18:  Re[Ex], Im[Ex], Re[Ey], Im[Ey], Re[Ez], Im[Ez]
        cols 19-24:  Re[dUx_ref], Im[dUx_ref], Re[dUy_ref], Im[dUy_ref], Re[dUz_ref], Im[dUz_ref]
    (1-based indexing shown here; see the table below for 0-based indices)

  - Klein & Howes (2015), arXiv:1503.00695:
      https://arxiv.org/pdf/1503.00695.pdf
    Used here for the standard (ω,γ) complex-frequency split and the sign interpretation:
      γ > 0 => unstable/growing mode; γ < 0 => stable/damped mode.

Helpful convention reference (plane wave notation):
  - MIT OCW plasma notes use wave quantities ∝ exp(i(k·x − ω t)) and take Re[...]:
      https://ocw.mit.edu/courses/22-611j-introduction-to-plasma-physics-i-fall-2003/75c75c5df51f8fba1b2a903fde937ddc_chap5.pdf

Assumptions (requested to note explicitly):
  - File encoding is UTF-8.
  - row_index is 0-based (first numeric data row is row_index=0).
  - time is dimensionless, i.e. (t * Omega_ref), consistent with PLUME output ω/Omega_ref.
  - Positions are in rho_ref units, so k_perp_rho*x and k_par_rho*z are dimensionless.

Implementation constraints:
  - Only uses: numpy, dataclasses, typing, pathlib
  - Optional visualization dependency: pyvista

-------------------------------------------------------------------------------
PLUME .modeN column map (for eigen=true), as a Markdown comment table:

| 0-based col | 1-based col | name (this script)                 |
|-------------|-------------|------------------------------------|
| 0           | 1           | k_perp_rho (= k_perp*rho_ref)      |
| 1           | 2           | k_par_rho  (= k_par*rho_ref)       |
| 2           | 3           | beta_ref_par (unused; in raw_row)  |
| 3           | 4           | vt_ref_par_over_c (unused; raw_row)|
| 4           | 5           | omega_r (= omega_r/Omega_ref)      |
| 5           | 6           | gamma   (= gamma/Omega_ref)        |
| 6           | 7           | Re[Bx]                             |
| 7           | 8           | Im[Bx]                             |
| 8           | 9           | Re[By]                             |
| 9           | 10          | Im[By]                             |
| 10          | 11          | Re[Bz]                             |
| 11          | 12          | Im[Bz]                             |
| 12          | 13          | Re[Ex]                             |
| 13          | 14          | Im[Ex]                             |
| 14          | 15          | Re[Ey]                             |
| 15          | 16          | Im[Ey]                             |
| 16          | 17          | Re[Ez]                             |
| 17          | 18          | Im[Ez]                             |
| 18          | 19          | Re[δU_ref,x]                       |
| 19          | 20          | Im[δU_ref,x]                       |
| 20          | 21          | Re[δU_ref,y]                       |
| 21          | 22          | Im[δU_ref,y]                       |
| 22          | 23          | Re[δU_ref,z]                       |
| 23          | 24          | Im[δU_ref,z]                       |

-------------------------------------------------------------------------------
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, List, Optional, Tuple, Union

import numpy as np

try:
    import pyvista as pv
except ImportError:  # pragma: no cover - exercised only when PyVista is absent.
    pv = None


@dataclass
class PlumeMode:
    """
    Minimal structured container for one PLUME eigenmode at one scan point.

    Required fields (per spec):
      - k_perp_rho, k_par_rho
      - omega_r, gamma
      - deltaB0: complex np.ndarray shape (3,) [Bx, By, Bz]
      - deltaU_ref0: complex np.ndarray shape (3,) [Ux, Uy, Uz]
      - raw_row: np.ndarray (all parsed floats from the selected row)

    Extra field:
      - unit_db_scale: cached normalization factor for unit_db=True
    """

    k_perp_rho: float
    k_par_rho: float
    omega_r: float
    gamma: float
    deltaB0: np.ndarray
    deltaU_ref0: np.ndarray
    raw_row: np.ndarray

    # Cached visualization scale factor (computed on first call to evaluate_deltaB with unit_db=True)
    unit_db_scale: Optional[float] = None


@dataclass
class DomainSampling:
    """
    Lightweight metadata describing the sampled visualization domain.
    """

    x_vals: np.ndarray
    z_vals: np.ndarray
    Lx: float
    Lz: float
    lambda_x: Optional[float]
    lambda_z: Optional[float]
    y: float


def load_mode(file_path: Union[str, Path], row_index: int) -> PlumeMode:
    """
    Read a PLUME .modeN file (ASCII, no header) and parse a single row by 0-based index.

    Parsing (minimum required):
      - k_perp_rho  = k_perp * rho_ref
      - k_par_rho   = k_par  * rho_ref
      - omega_r     = omega_r / Omega_ref  (dimensionless)
      - gamma       = gamma   / Omega_ref  (dimensionless)
      - deltaB0     = complex eigenvector components [Bx, By, Bz]
      - deltaU_ref0 = complex reference-species velocity components [Ux, Uy, Uz]

    Notes:
      - Uses line-by-line reading (memory-light).
      - Uses numpy.fromstring for fast whitespace parsing; sep=" " chosen explicitly.
      - Raises an informative error if eigen columns are missing (likely eigen=false run).

    See:
      https://github.com/kgklein/PLUME/blob/main/output.md
    """
    path = Path(file_path)

    if not isinstance(row_index, int):
        raise TypeError("row_index must be an int (0-based).")
    if row_index < 0:
        raise IndexError("row_index must be >= 0 (0-based indexing assumed).")
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    target_line: Optional[str] = None
    data_line_idx = -1  # counts only numeric data lines

    with path.open("r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                continue
            data_line_idx += 1
            if data_line_idx == row_index:
                target_line = s
                break

    if target_line is None:
        raise IndexError(
            f"Requested row_index={row_index}, but file contains only {data_line_idx + 1} data rows."
        )

    # IMPORTANT (NumPy behavior): pass sep explicitly so we are in TEXT mode, not binary mode.
    # See numpy docs: numpy.fromstring(..., sep=" ")
    row = np.fromstring(target_line, sep=" ", dtype=float)

    if row.size < 6:
        raise ValueError(
            f"Row {row_index} in {path} has only {row.size} columns; expected at least 6 "
            "for (k_perp_rho, k_par_rho, beta_ref_par, vt_ref_par/c, omega_r, gamma)."
        )

    # Required scan/solution values
    k_perp_rho = float(row[0])
    k_par_rho = float(row[1])
    omega_r = float(row[4])
    gamma = float(row[5])

    # Eigen δB columns require eigen=true in PLUME
    if row.size < 12:
        raise ValueError(
            f"Row {row_index} in {path} has only {row.size} columns; expected >=12 to parse "
            "δB eigenvector (Re/Im for Bx, By, Bz). "
            "Likely cause: PLUME run with eigen=false. "
            "See PLUME output.md for .modeN column ordering."
        )

    bx = row[6] + 1j * row[7]
    by = row[8] + 1j * row[9]
    bz = row[10] + 1j * row[11]
    deltaB0 = np.array([bx, by, bz], dtype=complex)

    if row.size < 24:
        raise ValueError(
            f"Row {row_index} in {path} has {row.size} columns; expected >=24 to parse "
            "reference-species velocity eigenvector δU_ref (Re/Im for Ux, Uy, Uz). "
            "See PLUME output.md for .modeN column ordering."
        )

    ux = row[18] + 1j * row[19]
    uy = row[20] + 1j * row[21]
    uz = row[22] + 1j * row[23]
    deltaU_ref0 = np.array([ux, uy, uz], dtype=complex)

    return PlumeMode(
        k_perp_rho=k_perp_rho,
        k_par_rho=k_par_rho,
        omega_r=omega_r,
        gamma=gamma,
        deltaB0=deltaB0,
        deltaU_ref0=deltaU_ref0,
        raw_row=row,
        unit_db_scale=None,
    )


def _evaluate_vector_field_raw(
    mode: PlumeMode,
    field0: np.ndarray,
    positions: np.ndarray,
    time: float,
    include_damping: bool,
) -> np.ndarray:
    """
    Internal helper: compute Re[field0 * exp(i(k·x - ω t))] without unit normalization.
    Returns real array shape (N,3).
    """
    pos = np.asarray(positions, dtype=float)
    field = np.asarray(field0, dtype=complex)
    if pos.ndim != 2 or pos.shape[1] != 3:
        raise ValueError(f"positions must have shape (N,3); got {pos.shape}.")
    if field.shape != (3,):
        raise ValueError(f"field0 must have shape (3,); got {field.shape}.")

    # Geometry: B0 = B0 z-hat; k = (k_perp, 0, k_par).
    # Positions are in rho_ref units, and k_{perp,par} are in (1/rho_ref) units via k*rho_ref,
    # so k·x is dimensionless:
    #
    #   k·x = (k_perp*rho_ref)*x + (k_par*rho_ref)*z
    #
    # with y uniform, ky=0.
    kdotx = mode.k_perp_rho * pos[:, 0] + mode.k_par_rho * pos[:, 2]

    # Phase convention (requested):
    #   exp(i(k·x - ω t))
    # with ω = ω_r (+ i*gamma if include_damping True).
    #
    # Under this convention, ω = ω_r + iγ implies exp(i(k·x - ω_r t)) * exp(γ t),
    # so γ<0 damps and γ>0 grows (consistent with sign usage in Klein & Howes 2015).
    omega = mode.omega_r + (1j * mode.gamma if include_damping else 0.0)
    phase = kdotx - omega * float(time)  # dimensionless
    osc = np.exp(1j * phase)

    # Broadcast complex eigenvector over all points, then take real part.
    vec = field[np.newaxis, :] * osc[:, np.newaxis]
    return np.real(vec)


def evaluate_deltaB(
    mode: PlumeMode,
    positions: np.ndarray,
    time: float,
    include_damping: bool = False,
    unit_db: bool = True,
) -> np.ndarray:
    """
    Compute the real magnetic fluctuation field δB(x,t) from a PLUME mode:

      δB(x,t) = Re[ δB0 * exp(i(k·x − ω t)) ]

    Inputs
    ------
    mode:
        PlumeMode dataclass from load_mode() or synthetic generator.
    positions:
        (N,3) array of positions in rho_ref units.
    time:
        scalar time in dimensionless units (t*Omega_ref), consistent with PLUME ω/Omega_ref.
    include_damping:
        If False (default), use ω = ω_r only (loopable animations).
        If True, use ω = ω_r + i*gamma (amplitude factor exp(gamma*t)).
    unit_db:
        If True (default), scale so max |δB| == 1 across the provided positions.
        The scale factor is computed ONCE at t=0 with include_damping=False and cached
        in mode.unit_db_scale. This keeps scaling consistent across frames.

    Returns
    -------
    deltaB:
        (N,3) real array.

    Notes on normalization
    ----------------------
    PLUME eigenfluctuations are normalized internally (per PLUME output.md, e.g., fields by Ex),
    so the absolute amplitude is not a physical amplitude. The unit_db option is purely a
    visualization helper that preserves polarization/phase relationships while normalizing
    the displayed magnitude.
    """
    if unit_db and mode.unit_db_scale is None:
        # Requested: compute normalization at t=0 without damping and cache it.
        db0 = _evaluate_vector_field_raw(
            mode,
            mode.deltaB0,
            positions,
            time=0.0,
            include_damping=False,
        )
        mags = np.linalg.norm(db0, axis=1)
        maxmag = float(np.max(mags)) if mags.size else 0.0

        if not np.isfinite(maxmag) or maxmag <= 0.0:
            # Degenerate case: identically zero field on this sampling.
            mode.unit_db_scale = 1.0
        else:
            mode.unit_db_scale = 1.0 / maxmag

    db = _evaluate_vector_field_raw(
        mode,
        mode.deltaB0,
        positions,
        time=time,
        include_damping=include_damping,
    )
    if unit_db:
        db = db * float(mode.unit_db_scale)
    return db


def evaluate_deltaU_ref(
    mode: PlumeMode,
    positions: np.ndarray,
    time: float,
    include_damping: bool = False,
) -> np.ndarray:
    """
    Compute the real reference-species velocity fluctuation field δU_ref(x,t):

      δU_ref(x,t) = Re[ δU_ref0 * exp(i(k·x − ω t)) ]

    The returned field is not independently normalized so its displayed magnitude can
    share the same glyph scaling policy as δB.
    """
    return _evaluate_vector_field_raw(
        mode,
        mode.deltaU_ref0,
        positions,
        time=time,
        include_damping=include_damping,
    )


def build_domain(
    mode: PlumeMode,
    n_wavelengths_x: float,
    n_wavelengths_z: float,
    n_lines: int,
    n_points_per_line: int,
    y: float = 0.0,
) -> Tuple[np.ndarray, List[slice], DomainSampling]:
    """
    Build a simple Cartesian sampling domain for the eigenmode.

    Domain geometry (requested):
      - Straight "field lines" parallel to z (consistent with B0 || z-hat)
      - Lines evenly spaced in x
      - y fixed (default 0.0)
      - Points sampled evenly along each line in z

    Wavelengths (requested):
      λx = 2π / k_perp_rho
      λz = 2π / k_par_rho

    Handling k=0 safely:
      - If k_perp_rho == 0, the wave is uniform in x; use a finite fallback Lx.
      - If k_par_rho  == 0, the wave is uniform in z; use a finite fallback Lz.
      - This avoids an infinite wavelength from producing an unusable plotting domain.

    Returns:
      positions: (N,3) with N = n_lines * n_points_per_line
      line_slices: list of slice objects marking each line in the flat positions array
      domain: lightweight metadata for plotting and diagnostics
    """
    if n_lines <= 0:
        raise ValueError("n_lines must be >= 1.")
    if n_points_per_line <= 0:
        raise ValueError("n_points_per_line must be >= 1.")
    if n_wavelengths_x < 0 or n_wavelengths_z < 0:
        raise ValueError("n_wavelengths_x and n_wavelengths_z must be >= 0.")

    kx = float(mode.k_perp_rho)
    kz = float(mode.k_par_rho)

    # λ = 2π/|k| (use |k| so wavelength is positive even if k is negative).
    # If k==0, wavelength is formally infinite, so pick a finite fallback domain.
    lambda_x: Optional[float]
    if np.isclose(kx, 0.0):
        lambda_x = None
        Lx = 1.0 * max(1.0, float(n_wavelengths_x))  # arbitrary but finite
    else:
        lambda_x = 2.0 * np.pi / abs(kx)
        Lx = float(n_wavelengths_x) * lambda_x

    lambda_z: Optional[float]
    if np.isclose(kz, 0.0):
        lambda_z = None
        Lz = 1.0 * max(1.0, float(n_wavelengths_z))  # arbitrary but finite
    else:
        lambda_z = 2.0 * np.pi / abs(kz)
        Lz = float(n_wavelengths_z) * lambda_z

    # Evenly spaced line anchors in x
    if n_lines == 1:
        x_vals = np.array([0.0], dtype=float)
    else:
        x_vals = np.linspace(0.0, Lx, n_lines, endpoint=False).astype(float)

    # Evenly spaced sample points along each line in z
    if n_points_per_line == 1:
        z_vals = np.array([0.0], dtype=float)
    else:
        z_vals = np.linspace(0.0, Lz, n_points_per_line, endpoint=False).astype(float)

    N = n_lines * n_points_per_line
    positions = np.empty((N, 3), dtype=float)
    line_slices: List[slice] = []

    idx = 0
    for x in x_vals:
        positions[idx : idx + n_points_per_line, 0] = x
        positions[idx : idx + n_points_per_line, 1] = float(y)
        positions[idx : idx + n_points_per_line, 2] = z_vals
        line_slices.append(slice(idx, idx + n_points_per_line))
        idx += n_points_per_line

    domain = DomainSampling(
        x_vals=x_vals,
        z_vals=z_vals,
        Lx=float(Lx),
        Lz=float(Lz),
        lambda_x=lambda_x,
        lambda_z=lambda_z,
        y=float(y),
    )

    return positions, line_slices, domain


def make_synthetic_circular_mode(
    k_perp_rho: float = 1.0,
    k_par_rho: float = 1.0,
    omega_r: float = 1.0,
    gamma: float = 0.0,
    handedness: str = "R",
) -> PlumeMode:
    """
    Synthetic test mode generator: circularly polarized δB in the x–y plane.

    δB0 is chosen as:
      - Right-handed (default): Bx = 1, By = +i, Bz = 0
      - Left-handed:            Bx = 1, By = -i, Bz = 0

    δU_ref0 is chosen to be nonzero and visibly distinct from δB0 so the overlay is
    easy to inspect in both static and animated test visualizations.

    This is purely for visualization/debugging and does not attempt to satisfy any
    particular plasma dispersion relation.
    """
    handedness = handedness.strip().upper()
    phase = 1j if handedness.startswith("R") else -1j

    deltaB0 = np.array([1.0 + 0j, phase * 1.0, 0.0 + 0j], dtype=complex)
    deltaU_ref0 = np.array([0.6 + 0j, 0.0 + 0j, phase * 0.35], dtype=complex)

    # Populate a minimal raw_row consistent with the 0-based table in the header.
    raw_row = np.full(24, np.nan, dtype=float)
    raw_row[0] = float(k_perp_rho)
    raw_row[1] = float(k_par_rho)
    raw_row[4] = float(omega_r)
    raw_row[5] = float(gamma)
    raw_row[6] = np.real(deltaB0[0])
    raw_row[7] = np.imag(deltaB0[0])
    raw_row[8] = np.real(deltaB0[1])
    raw_row[9] = np.imag(deltaB0[1])
    raw_row[10] = np.real(deltaB0[2])
    raw_row[11] = np.imag(deltaB0[2])
    raw_row[18] = np.real(deltaU_ref0[0])
    raw_row[19] = np.imag(deltaU_ref0[0])
    raw_row[20] = np.real(deltaU_ref0[1])
    raw_row[21] = np.imag(deltaU_ref0[1])
    raw_row[22] = np.real(deltaU_ref0[2])
    raw_row[23] = np.imag(deltaU_ref0[2])

    return PlumeMode(
        k_perp_rho=float(k_perp_rho),
        k_par_rho=float(k_par_rho),
        omega_r=float(omega_r),
        gamma=float(gamma),
        deltaB0=deltaB0,
        deltaU_ref0=deltaU_ref0,
        raw_row=raw_row,
        unit_db_scale=None,
    )


def _safe_period(omega_r: float) -> float:
    """
    Dimensionless wave period T = 2π/|omega_r|, with a safe fallback if omega_r ~ 0.
    """
    if not np.isfinite(omega_r) or np.isclose(omega_r, 0.0):
        return 1.0
    return 2.0 * np.pi / abs(float(omega_r))


def _env_flag(name: str, default: bool) -> bool:
    """
    Parse a boolean from an environment variable.
    """
    value = os.getenv(name)
    if value is None:
        return default
    return value.strip().lower() in {"1", "true", "yes", "on"}


def _env_int(name: str, default: int) -> int:
    """
    Parse an int from an environment variable.
    """
    value = os.getenv(name)
    return default if value is None else int(value)


def _env_float(name: str, default: float) -> float:
    """
    Parse a float from an environment variable.
    """
    value = os.getenv(name)
    return default if value is None else float(value)


def _env_path(name: str, default: Optional[Path]) -> Optional[Path]:
    """
    Parse a filesystem path from an environment variable.
    """
    value = os.getenv(name)
    if value is None or value.strip() == "":
        return default
    return Path(value).expanduser()


def _format_mode_value(value: float) -> str:
    """
    Format a mode scalar compactly for human-readable titles.
    """
    return f"{float(value):.6g}"


def _format_filename_value(value: float) -> str:
    """
    Format a mode scalar into a deterministic filename-safe token.
    """
    token = f"{float(value):.6g}"
    token = token.replace("-", "m")
    token = token.replace(".", "p")
    token = token.replace("+", "")
    return token


def _sanitize_filename_token(text: str) -> str:
    """
    Convert a user-facing label into a conservative filename component.
    """
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", text.strip())
    return cleaned.strip("_")


def _build_mode_metadata_title(mode: PlumeMode, title_keyword: str = "") -> str:
    """
    Build a compact scientific title string from the active mode metadata.
    """
    rho_ref_label = "ρᴿ"
    pieces: List[str] = []
    keyword = title_keyword.strip()
    if keyword:
        pieces.append(keyword)
    pieces.append(
        f"k⊥{rho_ref_label}={_format_mode_value(mode.k_perp_rho)}, "
        f"k∥{rho_ref_label}={_format_mode_value(mode.k_par_rho)}"
    )
    return " | ".join(pieces)


def _add_orientation_axes(plotter: Any) -> Any:
    """
    Add the small orientation-axis widget without changing its default labels.
    """
    return plotter.add_axes()


def _add_plot_bounds_with_rho_labels(plotter: Any) -> Any:
    """
    Add visible plot bounds/grid labels in rho_ref-normalized coordinates.
    """
    rho_ref_label = "rho_R"
    return plotter.show_grid(
        color="lightgray",
        xtitle=f"x/{rho_ref_label}",
        ytitle=f"y/{rho_ref_label}",
        ztitle=f"z/{rho_ref_label}",
    )


def _add_field_arrow_legend(
    plotter: Any,
    glyph_color_b: str,
    glyph_color_u: str,
    show_velocity: bool,
) -> Any:
    """
    Add a fixed legend box mapping arrow colors to the plotted fluctuation fields.
    """
    legend_entries = [["δB", glyph_color_b]]
    if show_velocity:
        legend_entries.append(["δUᴿ", glyph_color_u])
    return plotter.add_legend(
        labels=legend_entries,
        bcolor="white",
        border=True,
        face=None,
        loc="upper right",
        size=(0.16, 0.14),
    )


def _resolve_output_path_with_metadata(
    output_path: Union[str, Path],
    mode: PlumeMode,
    title_keyword: str = "",
) -> Path:
    """
    Append keyword and mode metadata to an export path while preserving its suffix.
    """
    path = Path(output_path)
    parts = [path.stem]

    keyword = _sanitize_filename_token(title_keyword)
    if keyword:
        parts.append(keyword)

    parts.append(f"kperp_{_format_filename_value(mode.k_perp_rho)}")
    parts.append(f"kpar_{_format_filename_value(mode.k_par_rho)}")

    resolved_name = "_".join(part for part in parts if part)
    return path.with_name(f"{resolved_name}{path.suffix}")


def _require_pyvista() -> Any:
    """
    Return the imported PyVista module or raise a clear error if unavailable.
    """
    if pv is None:
        raise ImportError(
            "PyVista is required for visualization but is not installed in this environment. "
            "Install it with 'pip install pyvista' or 'pip install pyvista imageio pillow'."
        )
    return pv


def _line_polydata_from_positions(positions: np.ndarray, line_slices: List[slice]) -> Any:
    """
    Build a PolyData containing polyline cells for the straight background field lines.
    """
    pv_mod = _require_pyvista()
    n_lines = len(line_slices)
    cells: List[int] = []

    for sl in line_slices:
        line_indices = np.arange(sl.start, sl.stop, dtype=np.int64)
        cells.extend([line_indices.size, *line_indices.tolist()])

    poly = pv_mod.PolyData(np.asarray(positions, dtype=float))
    poly.lines = np.asarray(cells, dtype=np.int64)
    poly.field_data["n_lines"] = np.array([n_lines], dtype=np.int64)
    return poly


def _subsample_indices(n_points: int, arrow_stride: int) -> np.ndarray:
    """
    Choose evenly spaced glyph sample indices from the full point set.
    """
    if n_points <= 0:
        return np.empty(0, dtype=int)
    if arrow_stride <= 0:
        raise ValueError("arrow_stride must be >= 1.")
    return np.arange(0, n_points, arrow_stride, dtype=int)


def _compute_glyph_factor(
    positions: np.ndarray,
    arrow_scale: float,
    domain: Optional[DomainSampling],
) -> float:
    """
    Convert a unit-normalized vector magnitude into a readable physical glyph length.
    """
    if arrow_scale <= 0.0:
        raise ValueError("arrow_scale must be > 0.")

    if domain is not None:
        base_length = max(domain.Lx, domain.Lz, 1.0)
    elif positions.size:
        extents = np.ptp(positions, axis=0)
        base_length = max(float(np.max(extents)), 1.0)
    else:
        base_length = 1.0

    return float(arrow_scale) * base_length


def _default_camera_position(domain: Optional[DomainSampling]) -> Tuple[Tuple[float, float, float], Tuple[float, float, float], Tuple[float, float, float]]:
    """
    Deterministic oblique camera that shows x and z structure without any animation rotation.
    """
    if domain is None:
        center = (0.0, 0.0, 0.0)
        span = 1.0
    else:
        center = (0.5 * domain.Lx, domain.y, 0.5 * domain.Lz)
        span = max(domain.Lx, domain.Lz, 1.0)

    position = (center[0] + 1.8 * span, center[1] + 1.1 * span, center[2] + 1.6 * span)
    viewup = (0.0, 1.0, 0.0)
    return position, center, viewup


def _build_glyph_mesh(
    positions: np.ndarray,
    vectors: np.ndarray,
    sample_idx: np.ndarray,
    glyph_factor: float,
) -> Any:
    """
    Build a glyph mesh for a sampled vector field.
    """
    pv_mod = _require_pyvista()
    glyph_points = pv_mod.PolyData(np.asarray(positions, dtype=float)[sample_idx])
    sampled_vectors = np.asarray(vectors, dtype=float)[sample_idx]
    glyph_points["vectors"] = sampled_vectors
    glyph_points["vector_mag"] = np.linalg.norm(sampled_vectors, axis=1)
    return glyph_points.glyph(
        orient="vectors",
        scale="vector_mag",
        factor=glyph_factor,
        geom=pv_mod.Arrow(),
    )


def create_static_plot(
    mode: PlumeMode,
    positions: np.ndarray,
    line_slices: List[slice],
    deltaB_vectors: np.ndarray,
    deltaU_ref_vectors: Optional[np.ndarray] = None,
    domain: Optional[DomainSampling] = None,
    arrow_stride: int = 8,
    arrow_scale: float = 0.12,
    line_color: str = "lightgray",
    glyph_color_b: str = "royalblue",
    glyph_color_u: str = "gold",
    line_width: float = 1.0,
    show_velocity: bool = True,
    show: bool = True,
    window_size: Tuple[int, int] = (1200, 900),
    title_keyword: str = "",
) -> Any:
    """
    Create a static PyVista 3D scene showing background field lines and δB glyph arrows.

    The arrow length scale is global and fixed for the whole domain, so all glyphs share
    the same magnitude-to-length mapping.
    """
    pv_mod = _require_pyvista()
    pos = np.asarray(positions, dtype=float)
    db = np.asarray(deltaB_vectors, dtype=float)
    du = None if deltaU_ref_vectors is None else np.asarray(deltaU_ref_vectors, dtype=float)

    if pos.ndim != 2 or pos.shape[1] != 3:
        raise ValueError(f"positions must have shape (N,3); got {pos.shape}.")
    if db.shape != pos.shape:
        raise ValueError(f"deltaB_vectors must have shape {pos.shape}; got {db.shape}.")
    if du is not None and du.shape != pos.shape:
        raise ValueError(f"deltaU_ref_vectors must have shape {pos.shape}; got {du.shape}.")

    sample_idx = _subsample_indices(pos.shape[0], arrow_stride)
    glyph_factor = _compute_glyph_factor(pos, arrow_scale=arrow_scale, domain=domain)

    line_mesh = _line_polydata_from_positions(pos, line_slices)
    b_glyphs = _build_glyph_mesh(pos, db, sample_idx, glyph_factor)
    u_glyphs = None
    if show_velocity and du is not None:
        u_glyphs = _build_glyph_mesh(pos, du, sample_idx, glyph_factor)

    plotter = pv_mod.Plotter(window_size=window_size)
    plotter.add_mesh(line_mesh, color=line_color, line_width=float(line_width), render_lines_as_tubes=False)
    plotter.add_mesh(b_glyphs, color=glyph_color_b)
    if u_glyphs is not None:
        plotter.add_mesh(u_glyphs, color=glyph_color_u)
    _add_field_arrow_legend(
        plotter,
        glyph_color_b=glyph_color_b,
        glyph_color_u=glyph_color_u,
        show_velocity=u_glyphs is not None,
    )
    _add_orientation_axes(plotter)
    plotter.set_background("white")
    plotter.camera_position = _default_camera_position(domain)
    _add_plot_bounds_with_rho_labels(plotter)

    mode_metadata = _build_mode_metadata_title(mode, title_keyword=title_keyword)
    title = f"δB + δUᴿ static view | {mode_metadata}"
    plotter.add_title(title, font_size=12)

    if show:
        plotter.show()

    return plotter


def create_animation(
    mode: PlumeMode,
    positions: np.ndarray,
    line_slices: List[slice],
    output_gif_path: Union[str, Path],
    domain: Optional[DomainSampling] = None,
    frames_per_period: int = 16,
    arrow_stride: int = 8,
    arrow_scale: float = 0.12,
    include_damping: bool = False,
    line_color: str = "lightgray",
    glyph_color_b: str = "royalblue",
    glyph_color_u: str = "gold",
    line_width: float = 1.0,
    window_size: Tuple[int, int] = (1200, 900),
    show_progress: bool = True,
    show_velocity: bool = True,
    title_keyword: str = "",
) -> Path:
    """
    Generate a loopable PyVista GIF animation of δB over one wave period.

    Frames span [0, T) with T = 2π/|ωr|, so the wrap from the last frame back to the first
    is smooth without duplicating the initial state.
    """
    pv_mod = _require_pyvista()
    pos = np.asarray(positions, dtype=float)
    if pos.ndim != 2 or pos.shape[1] != 3:
        raise ValueError(f"positions must have shape (N,3); got {pos.shape}.")
    if frames_per_period < 2:
        raise ValueError("frames_per_period must be >= 2 for a loopable animation.")

    output_path = _resolve_output_path_with_metadata(
        output_gif_path,
        mode=mode,
        title_keyword=title_keyword,
    )
    sample_idx = _subsample_indices(pos.shape[0], arrow_stride)
    glyph_factor = _compute_glyph_factor(pos, arrow_scale=arrow_scale, domain=domain)
    line_mesh = _line_polydata_from_positions(pos, line_slices)

    period = _safe_period(mode.omega_r)
    times = np.linspace(0.0, period, frames_per_period, endpoint=False)
    total_frames = times.size

    if show_progress:
        print(
            f"Starting GIF render: {output_path} "
            f"({total_frames} frames, include_damping={include_damping}, "
            f"show_velocity={show_velocity})"
        )

    db0 = evaluate_deltaB(mode, pos, time=0.0, include_damping=include_damping, unit_db=True)
    b_glyphs = _build_glyph_mesh(pos, db0, sample_idx, glyph_factor)
    u_glyphs = None
    if show_velocity:
        du0 = evaluate_deltaU_ref(mode, pos, time=0.0, include_damping=include_damping)
        u_glyphs = _build_glyph_mesh(pos, du0, sample_idx, glyph_factor)

    plotter = pv_mod.Plotter(off_screen=True, window_size=window_size)
    plotter.open_gif(str(output_path))
    plotter.add_mesh(line_mesh, color=line_color, line_width=float(line_width), render_lines_as_tubes=False)
    b_actor = plotter.add_mesh(b_glyphs, color=glyph_color_b)
    u_actor = None
    if u_glyphs is not None:
        u_actor = plotter.add_mesh(u_glyphs, color=glyph_color_u)
    _add_field_arrow_legend(
        plotter,
        glyph_color_b=glyph_color_b,
        glyph_color_u=glyph_color_u,
        show_velocity=u_glyphs is not None,
    )
    _add_orientation_axes(plotter)
    plotter.set_background("white")
    plotter.camera_position = _default_camera_position(domain)
    _add_plot_bounds_with_rho_labels(plotter)

    mode_metadata = _build_mode_metadata_title(mode, title_keyword=title_keyword)
    title = f"δB + δUᴿ animation | {mode_metadata}"
    plotter.add_title(title, font_size=12)

    for frame_idx, time in enumerate(times, start=1):
        if show_progress:
            print(f"Rendering frame {frame_idx}/{total_frames} (t = {float(time):.6f})")

        db = evaluate_deltaB(mode, pos, time=float(time), include_damping=include_damping, unit_db=True)
        b_frame_glyphs = _build_glyph_mesh(pos, db, sample_idx, glyph_factor)

        plotter.remove_actor(b_actor, reset_camera=False)
        b_actor = plotter.add_mesh(b_frame_glyphs, color=glyph_color_b)

        if show_velocity:
            du = evaluate_deltaU_ref(mode, pos, time=float(time), include_damping=include_damping)
            u_frame_glyphs = _build_glyph_mesh(pos, du, sample_idx, glyph_factor)
            if u_actor is not None:
                plotter.remove_actor(u_actor, reset_camera=False)
            u_actor = plotter.add_mesh(u_frame_glyphs, color=glyph_color_u)

        plotter.write_frame()

    plotter.close()
    if show_progress:
        print(f"Finished GIF render: {output_path}")
    return output_path


if __name__ == "__main__":
    # -------------------------------------------------------------------------
    # User-editable demo configuration.
    # These defaults can also be overridden from the shell with environment
    # variables such as SWIFT_CREATE_STATIC_PLOT=1 or SWIFT_OUTPUT_GIF=foo.gif.
    # -------------------------------------------------------------------------

    # Example PLUME file (edit to your local path). If missing, script uses synthetic mode.
    MODE_FILE: Optional[Path] = None
    # MODE_FILE = Path("data/example/map_par_kpar_1_10000.mode1")
    MODE_FILE = _env_path("SWIFT_MODE_FILE", MODE_FILE)

    ROW_INDEX = _env_int("SWIFT_ROW_INDEX", 0)

    # Domain controls
    N_WAVELENGTHS_X = _env_float("SWIFT_N_WAVELENGTHS_X", 1.0)
    N_WAVELENGTHS_Z = _env_float("SWIFT_N_WAVELENGTHS_Z", 1.0)
    N_LINES = _env_int("SWIFT_N_LINES", 6)
    N_POINTS_PER_LINE = _env_int("SWIFT_N_POINTS_PER_LINE", 64)
    Y0 = _env_float("SWIFT_Y0", 0.0)

    # Visualization controls
    CREATE_STATIC_PLOT = _env_flag("SWIFT_CREATE_STATIC_PLOT", False)
    CREATE_GIF = _env_flag("SWIFT_CREATE_GIF", False)
    OUTPUT_GIF = _env_path("SWIFT_OUTPUT_GIF", Path("deltaB_mode.gif"))
    ARROW_STRIDE = _env_int("SWIFT_ARROW_STRIDE", 8)
    ARROW_SCALE = _env_float("SWIFT_ARROW_SCALE", 0.12)
    FRAMES_PER_PERIOD = _env_int("SWIFT_FRAMES_PER_PERIOD", 16)
    SHOW_PROGRESS = _env_flag("SWIFT_SHOW_PROGRESS", True)
    SHOW_VELOCITY = _env_flag("SWIFT_SHOW_VELOCITY", True)
    TITLE_KEYWORD = os.getenv("SWIFT_TITLE_KEYWORD", "").strip()

    # -------------------------------------------------------------------------
    # Load a mode (PLUME or synthetic)
    # -------------------------------------------------------------------------
    if MODE_FILE is not None and MODE_FILE.exists():
        mode = load_mode(MODE_FILE, ROW_INDEX)
        print(f"Loaded PLUME mode row {ROW_INDEX} from: {MODE_FILE}")
    else:
        mode = make_synthetic_circular_mode(
            k_perp_rho=1.0,
            k_par_rho=0.5,
            omega_r=1.0,
            gamma=0.0,
            handedness="R",
        )
        print("MODE_FILE not provided or not found; using synthetic circularly polarized test mode.")

    print(f"k_perp_rho={mode.k_perp_rho:.6g}, k_par_rho={mode.k_par_rho:.6g}")
    print(f"omega_r={mode.omega_r:.6g}, gamma={mode.gamma:.6g}")
    print(f"deltaB0 (complex) = {mode.deltaB0}")
    print(f"deltaU_ref0 (complex) = {mode.deltaU_ref0}")

    # -------------------------------------------------------------------------
    # Build a sampling domain and evaluate δB at two times
    # -------------------------------------------------------------------------
    positions, line_slices, domain = build_domain(
        mode=mode,
        n_wavelengths_x=N_WAVELENGTHS_X,
        n_wavelengths_z=N_WAVELENGTHS_Z,
        n_lines=N_LINES,
        n_points_per_line=N_POINTS_PER_LINE,
        y=Y0,
    )

    print(f"positions shape = {positions.shape} (N = {positions.shape[0]})")
    print(f"number of field lines = {len(line_slices)}; points per line = {N_POINTS_PER_LINE}")
    print(f"Lx = {domain.Lx:.6g}, Lz = {domain.Lz:.6g}")
    print(f"lambda_x = {domain.lambda_x}, lambda_z = {domain.lambda_z}")

    # Evaluate at t=0
    db_t0 = evaluate_deltaB(mode, positions, time=0.0, include_damping=False, unit_db=True)
    du_t0 = evaluate_deltaU_ref(mode, positions, time=0.0, include_damping=False)
    maxmag_t0 = np.max(np.linalg.norm(db_t0, axis=1))
    maxu_t0 = np.max(np.linalg.norm(du_t0, axis=1))
    print(f"deltaB(t=0) shape = {db_t0.shape}")
    print(f"deltaU_ref(t=0) shape = {du_t0.shape}")
    print(f"unit_db_scale (cached) = {mode.unit_db_scale}")
    print(f"max |deltaB(t=0)| across positions = {maxmag_t0:.6g} (should be 1.0 if nonzero)")
    print(f"max |deltaU_ref(t=0)| across positions = {maxu_t0:.6g}")

    # Evaluate one frame later (dt = period / frames)
    period = _safe_period(mode.omega_r)
    dt = period / float(FRAMES_PER_PERIOD)
    db_t1 = evaluate_deltaB(mode, positions, time=dt, include_damping=False, unit_db=True)
    du_t1 = evaluate_deltaU_ref(mode, positions, time=dt, include_damping=False)
    maxmag_t1 = np.max(np.linalg.norm(db_t1, axis=1))
    maxu_t1 = np.max(np.linalg.norm(du_t1, axis=1))
    print(f"deltaB(t=dt) shape = {db_t1.shape}, dt = {dt:.6g}")
    print(f"deltaU_ref(t=dt) shape = {du_t1.shape}, dt = {dt:.6g}")
    print(f"max |deltaB(t=dt)| across positions = {maxmag_t1:.6g} (should remain ~1.0)")
    print(f"max |deltaU_ref(t=dt)| across positions = {maxu_t1:.6g}")

    # -------------------------------------------------------------------------
    # Optional PyVista visualization
    # -------------------------------------------------------------------------
    if CREATE_STATIC_PLOT:
        try:
            create_static_plot(
                mode=mode,
                positions=positions,
                line_slices=line_slices,
                deltaB_vectors=db_t0,
                deltaU_ref_vectors=du_t0,
                domain=domain,
                arrow_stride=ARROW_STRIDE,
                arrow_scale=ARROW_SCALE,
                show_velocity=SHOW_VELOCITY,
                show=True,
                title_keyword=TITLE_KEYWORD,
            )
        except Exception as e:
            print(f"Static PyVista plot requested but failed: {e}")

    if CREATE_GIF:
        try:
            gif_path = create_animation(
                mode=mode,
                positions=positions,
                line_slices=line_slices,
                output_gif_path=OUTPUT_GIF,
                domain=domain,
                frames_per_period=FRAMES_PER_PERIOD,
                arrow_stride=ARROW_STRIDE,
                arrow_scale=ARROW_SCALE,
                include_damping=False,
                show_progress=SHOW_PROGRESS,
                show_velocity=SHOW_VELOCITY,
                title_keyword=TITLE_KEYWORD,
            )
            print(f"Wrote loopable GIF to: {gif_path}")
        except Exception as e:
            print(f"GIF animation requested but failed: {e}")
