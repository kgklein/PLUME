#!/usr/bin/env python3
"""
plume_wave_viz.py

Single-file utilities to:
  - Parse a single eigenmode row from a PLUME ASCII .modeN file (no header).
  - Reconstruct the real-space magnetic fluctuation field \delta B(x,t) from the
    complex eigenvector and complex frequency.
  - Build a simple Cartesian sampling domain aligned with B0 = B0 z-hat.

No PyVista plotting is included in this step (per current scope).
This script is intended to be readable and easy to modify.

Primary format & convention references (as requested):
  - PLUME output column ordering (.modeN files):
      https://github.com/kgklein/PLUME/blob/main/output.md
    In particular, for .modeN from om_scan:
      cols 1-6: k_perp*rho_ref, k_par*rho_ref, beta_ref_par, vt_ref_par/c, omega_r/Omega_ref, gamma/Omega_ref
      if eigen=true, cols 7-12: Re[Bx], Im[Bx], Re[By], Im[By], Re[Bz], Im[Bz]
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
  - Only uses: numpy, dataclasses, typing, pathlib (plus optional matplotlib in a debug block).

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

-------------------------------------------------------------------------------
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple, Union

import numpy as np


@dataclass
class PlumeMode:
    """
    Minimal structured container for one PLUME eigenmode at one scan point.

    Required fields (per spec):
      - k_perp_rho, k_par_rho
      - omega_r, gamma
      - deltaB0: complex np.ndarray shape (3,) [Bx, By, Bz]
      - raw_row: np.ndarray (all parsed floats from the selected row)

    Extra field:
      - unit_db_scale: cached normalization factor for unit_db=True
    """

    k_perp_rho: float
    k_par_rho: float
    omega_r: float
    gamma: float
    deltaB0: np.ndarray
    raw_row: np.ndarray

    # Cached visualization scale factor (computed on first call to evaluate_deltaB with unit_db=True)
    unit_db_scale: Optional[float] = None


def load_mode(file_path: Union[str, Path], row_index: int) -> PlumeMode:
    """
    Read a PLUME .modeN file (ASCII, no header) and parse a single row by 0-based index.

    Parsing (minimum required):
      - k_perp_rho  = k_perp * rho_ref
      - k_par_rho   = k_par  * rho_ref
      - omega_r     = omega_r / Omega_ref  (dimensionless)
      - gamma       = gamma   / Omega_ref  (dimensionless)
      - deltaB0     = complex eigenvector components [Bx, By, Bz]

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

    return PlumeMode(
        k_perp_rho=k_perp_rho,
        k_par_rho=k_par_rho,
        omega_r=omega_r,
        gamma=gamma,
        deltaB0=deltaB0,
        raw_row=row,
        unit_db_scale=None,
    )


def _deltaB_raw(
    mode: PlumeMode,
    positions: np.ndarray,
    time: float,
    include_damping: bool,
) -> np.ndarray:
    """
    Internal helper: compute δB(x,t) without unit normalization.
    Returns real array shape (N,3).
    """
    pos = np.asarray(positions, dtype=float)
    if pos.ndim != 2 or pos.shape[1] != 3:
        raise ValueError(f"positions must have shape (N,3); got {pos.shape}.")

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
    vec = mode.deltaB0[np.newaxis, :] * osc[:, np.newaxis]
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
        db0 = _deltaB_raw(mode, positions, time=0.0, include_damping=False)
        mags = np.linalg.norm(db0, axis=1)
        maxmag = float(np.max(mags)) if mags.size else 0.0

        if not np.isfinite(maxmag) or maxmag <= 0.0:
            # Degenerate case: identically zero field on this sampling.
            mode.unit_db_scale = 1.0
        else:
            mode.unit_db_scale = 1.0 / maxmag

    db = _deltaB_raw(mode, positions, time=time, include_damping=include_damping)
    if unit_db:
        db = db * float(mode.unit_db_scale)
    return db


def build_domain(
    n_wavelengths_x: float,
    n_wavelengths_z: float,
    n_lines: int,
    n_points_per_line: int,
    k_perp_rho: float,
    k_par_rho: float,
    y: float = 0.0,
) -> Tuple[np.ndarray, List[slice]]:
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
      - If k_perp_rho == 0, the wave is uniform in x; choose an arbitrary finite Lx.
      - If k_par_rho  == 0, the wave is uniform in z; choose an arbitrary finite Lz.

    Returns:
      positions: (N,3) with N = n_lines * n_points_per_line
      line_slices: list of slice objects marking each line in the flat positions array
    """
    if n_lines <= 0:
        raise ValueError("n_lines must be >= 1.")
    if n_points_per_line <= 0:
        raise ValueError("n_points_per_line must be >= 1.")
    if n_wavelengths_x < 0 or n_wavelengths_z < 0:
        raise ValueError("n_wavelengths_x and n_wavelengths_z must be >= 0.")

    kx = float(k_perp_rho)
    kz = float(k_par_rho)

    # λ = 2π/|k| (use |k| so wavelength is positive even if k is negative).
    # If k==0, wavelength is effectively infinite, so pick a finite domain size.
    if np.isclose(kx, 0.0):
        Lx = 1.0 * max(1.0, float(n_wavelengths_x))  # arbitrary but finite
    else:
        lambda_x = 2.0 * np.pi / abs(kx)
        Lx = float(n_wavelengths_x) * lambda_x

    if np.isclose(kz, 0.0):
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

    return positions, line_slices


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

    This is purely for visualization/debugging and does not attempt to satisfy any
    particular plasma dispersion relation.
    """
    handedness = handedness.strip().upper()
    phase = 1j if handedness.startswith("R") else -1j

    deltaB0 = np.array([1.0 + 0j, phase * 1.0, 0.0 + 0j], dtype=complex)

    # Populate a minimal raw_row consistent with the 0-based table in the header.
    raw_row = np.full(12, np.nan, dtype=float)
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

    return PlumeMode(
        k_perp_rho=float(k_perp_rho),
        k_par_rho=float(k_par_rho),
        omega_r=float(omega_r),
        gamma=float(gamma),
        deltaB0=deltaB0,
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


if __name__ == "__main__":
    # -------------------------------------------------------------------------
    # User-editable demo configuration (no CLI in this step by design).
    # -------------------------------------------------------------------------

    # Example PLUME file (edit to your local path). If missing, script uses synthetic mode.
    MODE_FILE: Optional[Path] = None
    # MODE_FILE = Path("data/example/map_par_kpar_1_10000.mode1")

    ROW_INDEX = 0

    # Domain controls
    N_WAVELENGTHS_X = 1.0
    N_WAVELENGTHS_Z = 1.0
    N_LINES = 6
    N_POINTS_PER_LINE = 64
    Y0 = 0.0

    # Animation-relevant parameter used only for dt demonstration (no plotting here)
    FRAMES_PER_PERIOD = 16

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

    # -------------------------------------------------------------------------
    # Build a sampling domain and evaluate δB at two times
    # -------------------------------------------------------------------------
    positions, line_slices = build_domain(
        n_wavelengths_x=N_WAVELENGTHS_X,
        n_wavelengths_z=N_WAVELENGTHS_Z,
        n_lines=N_LINES,
        n_points_per_line=N_POINTS_PER_LINE,
        k_perp_rho=mode.k_perp_rho,
        k_par_rho=mode.k_par_rho,
        y=Y0,
    )

    print(f"positions shape = {positions.shape} (N = {positions.shape[0]})")
    print(f"number of field lines = {len(line_slices)}; points per line = {N_POINTS_PER_LINE}")

    # Evaluate at t=0
    db_t0 = evaluate_deltaB(mode, positions, time=0.0, include_damping=False, unit_db=True)
    maxmag_t0 = np.max(np.linalg.norm(db_t0, axis=1))
    print(f"deltaB(t=0) shape = {db_t0.shape}")
    print(f"unit_db_scale (cached) = {mode.unit_db_scale}")
    print(f"max |deltaB(t=0)| across positions = {maxmag_t0:.6g} (should be 1.0 if nonzero)")

    # Evaluate one frame later (dt = period / frames)
    period = _safe_period(mode.omega_r)
    dt = period / float(FRAMES_PER_PERIOD)
    db_t1 = evaluate_deltaB(mode, positions, time=dt, include_damping=False, unit_db=True)
    maxmag_t1 = np.max(np.linalg.norm(db_t1, axis=1))
    print(f"deltaB(t=dt) shape = {db_t1.shape}, dt = {dt:.6g}")
    print(f"max |deltaB(t=dt)| across positions = {maxmag_t1:.6g} (should remain ~1.0)")

    # -------------------------------------------------------------------------
    # Optional quick matplotlib debug (disabled by default).
    # This plots arrows in x–z using (Bx, Bz) at y=constant.
    # -------------------------------------------------------------------------
    USE_MATPLOTLIB_DEBUG = False
    if USE_MATPLOTLIB_DEBUG:
        try:
            import matplotlib.pyplot as plt  # optional

            x = positions[:, 0]
            z = positions[:, 2]
            bx = db_t0[:, 0]
            bz = db_t0[:, 2]

            plt.figure()
            plt.quiver(x, z, bx, bz, angles="xy", scale_units="xy", scale=1.0)
            plt.xlabel("x / rho_ref")
            plt.ylabel("z / rho_ref")
            plt.title("δB quiver (Bx,Bz) at t=0 (unit_db normalized)")
            plt.axis("equal")
            plt.show()
        except Exception as e:
            print(f"Matplotlib debug requested but failed: {e}")
