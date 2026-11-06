<a id="jet-plume-historical-field-particle-correlations"></a>

# JET-PLUME: historical field-particle correlations

[Manual home](./) · [Input reference](input.md#fpc-fpc) · [Output readers](output.md#fpc-and-debug-families)

JET-PLUME (Judging Energy Transfer in a Plasma in a Linear Uniform Magnetized
Environment) reconstructs linear velocity-space distribution perturbations and
field-particle correlation (FPC) signatures. It is an **all-kinetic historical
pathway**, separate from corrected species Ps and wave energy. Cold/mixed FPC is
not implemented, including warm populations in a mixture. New density/velocity
moments and corrected heating do not enable that future functionality.

Original JET-PLUME contributions are credited to Collin Brown, Gregory Howes,
Kristopher Klein and Jason TenBarge. The [PLUME paper and existing scientific
references](./#authors-scientific-references-and-citation) provide context.

<a id="workflows-and-coordinates"></a>

## Workflows and coordinates

`option=6` invokes `compute_fpc_gyro`: distributions/correlations are reduced over
velocity gyrophase to a (vperp,vpar) grid. `option=7` invokes `compute_fpc_cart`:
Cartesian velocity-space reconstruction, three projected planes, and 3D/debug
outputs. Earlier prose reversing these routines was incorrect. In both cases
use_map selects map initialization versus supplied guesses refined by the solver;
there is not an independent mixed/cold FPC root solver.

The background field is z, k lies in x-z, x is perpendicular in the wave/field
plane, y is its second perpendicular direction, and z is parallel. Grid velocities
use the historical species **parallel thermal-speed** coordinates, not the
species-moment output normalization. FPC grid construction/header values still
use the derived rho_ref solver convention for d_ref inputs.

Add `&fpc` to an otherwise valid all-kinetic input. This is a parameter fragment,
not a standalone runnable plasma:

```fortran
&fpc
vperpmin=0, vperpmax=.5, vparmin=-.5, vparmax=.5
vxmin=-.5, vxmax=.5, vymin=-.5, vymax=.5, vzmin=-.5, vzmax=.5
delv=.5, elecdircontribution=0
/
```

Those bounds/spacing are exercised by the existing FPC normalization test:

```bash
PYTHONPATH=tests MPLCONFIGDIR=/tmp/plume-mpl python3 -m unittest \
  test_normalization.NormalizationTests.test_fpc_diagnostics -v
```

A curated standalone public FPC lesson is still outstanding. The existing
[example notebook](https://github.com/kgklein/PLUME/blob/main/examplelinfpc.ipynb)
shows the optional wrapper's historical workflow; review mode/schema assumptions
before reuse. This stage does not create or rewrite input examples.

Use positive delv and integer-compatible bounds. The writer derives integer grid
indices; noninteger-multiple bounds warn and may cause off-by-one reader behavior.
Gyro bounds matter for option6; Cartesian bounds for option7. delv is shared across
directions. elecdircontribution is a debug forcing selector: 1 keeps Ex, 2 Ey,
3 Ez, and 0/other leaves forcing unchanged. It is not a heating-channel selector.

<a id="meaning-and-normalization-caveats"></a>

## Meaning and normalization caveats

The distribution files contain real/imaginary parts of the internally normalized
linear Fourier coefficient f_s1, not an independent physical measured distribution.
For Cartesian internal thermal-coordinate velocity vhat_i and distribution
coefficient fhat_1, the literal implemented weight is

\[
 \widehat C_i=-\frac{Q_s}{2}\widehat v_i^2
 \operatorname{Re}\!\left[E_i^*
                 \frac{\partial\widehat f_1}{\partial\widehat v_i}\right].
\]

Here Q_s is **the code's input charge ratio** q_R/q_s, E_i is the Ex=1 field
component, and vhat/fhat are the FPC routine's internal thermal variables. This
states the implemented nondimensional weight, not a new dimensional-power
normalization. Projection integrates the omitted coordinate with the grid
weights; gyro reduction includes the cylindrical Jacobian/gyrophase weights
and a special near-zero radial cell treatment. There are existing finite-grid,
harmonic and normalization limitations. Do not equate these raw weights or their
debug integrals with corrected Ps or entropy production without a separate
normalization/closure analysis. A universal dimensional conversion for these
historical FPC outputs is not newly validated by this manual.

Headers labeled omega.r/omega.i literally contain real/imaginary parts of
w sqrt(beta_R)/Kpar, where w=omega/Omega_R and Kpar=kpar rho_R in solver
coordinates. Since rho_R/d_R=sqrt(beta_R alpha_R), this equals
\((w/\kappa_\parallel)/\sqrt{\alpha_R}\), **not** the physical phase speed in
v_A,R units when alpha_R differs from one. Old stdout/plot labels claiming an
unqualified Alfvén-normalized phase speed need that caveat. A future FPC audit
must address the full header/grid/normalization contract before broader use.

FPC reconstruction has fixed internal harmonic range nbesmax=15; adaptive
susceptibility Bessels do not change it. Numerical distribution moments/debug
integrals are not the validated analytic current/continuity species moments.
Their quadrature can be inaccurate. No mixed FPC or new collisional accuracy
claim is made here.

<a id="files-and-array-orientation"></a>

## Files and array orientation

One indexed set is produced per population and root, with species in the original
input order, generally `.specieSS.modeMM`. Old readers use index-based electron
labels and offset recovery; those heuristics are not arbitrary-mixture validation.

<a id="gyro-files-option6"></a>

### Gyro files: option6

| Quantity suffix | Contents |
|---|---|
| .cpar | Parallel-field correlation |
| .cperp | Perpendicular-field correlation |
| .df1gyro.real | Real reduced perturbed distribution coefficient |
| .df1gyro.imag | Imaginary reduced perturbed distribution coefficient |

Each has eight physical/thermal header labels/values, five grid labels/values,
a separator, then a matrix. vperp advances by row, vpar by column; the first
entry corresponds to the minimum grid values. The matrix is a gyro-reduced
velocity-space object, not a projected Cartesian map mislabeled as gyro.

<a id="cartesian-files-option7"></a>

### Cartesian files: option7

| Quantity suffix | Contents |
|---|---|
| .cparcart | Ez correlation, projected into three planes |
| .cperp1 | Ex correlation, projected into three planes |
| .cperp2 | Ey correlation, projected into three planes |
| .dfs.real, .dfs.imag | Real/imaginary distribution projected into those planes |
| .dfs3D.real, .dfs3D.imag | Real/imaginary full 3D distribution blocks |

Projection files contain (x,y), (x,z), (z,y) array layouts separated by `---`.
The first has rows x / columns y; second rows x / columns z; third rows z /
columns y. This last historical orientation assumes compatible y/z grids;
it is not a reason to transpose every loaded diagnostic blindly. The old
`.df1.real/.imag` naming in prose does not match the current `.dfs` writer.

Three-dimensional blocks use reader storage (z,y,x), one writer line per (x,y)
containing z values. They have a separate block contract and must not be passed
to the projected-file parser. Source/header-specific debug `.eigen`/`.mom` files
are described on [Output](output.md#fpc-and-debug-families); they have no ordinary
scan species-parameter suffix and are unsuitable as a general scientific reader.

Some FPC **writers** already zero extremely small values and saturate large ones
(typical thresholds 9.999e-99 and 9.999e99). Reader corrections cannot reconstruct
values lost by that historical formatting. A zero in such a file is not proof of
an exactly vanishing physical correlation; inspect its source/precision regime.
No clipping thresholds or numerical FPC calculation are changed by this manual.

<a id="reading-and-missing-values"></a>

## Reading and missing values

Use `loadlinfpccepar`, `loadlinfpcceperp`, `loadlinfpcgyro_dist` for gyro;
`loadlinfpccart`/`loadlinfpccart_dist` for projections;
`loadlinfpc3d`/`loadlinfpc3d_dist` for full Cartesian distribution blocks.
Readers preserve NaNs, infinities, genuine finite zeros, signs and shape by
default. `visualization_nan_to_zero=True` explicitly replaces **only NaNs** for
plots; it is not a physical value or analysis correction. The 3D reader no
longer clips representable finite samples itself, though old writers may already
have done so. Distinguish a missing diagnostic from a physical calculated zero.
