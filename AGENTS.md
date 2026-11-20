
# PLUME — Repository Instructions

## Project overview

PLUME is a Fortran plasma dispersion relation solver for linear waves and instabilities in magnetized plasmas. It supports kinetic and cold-plasma species, parallel drifts, and oblique propagation.

The primary development priorities are scientific correctness, numerical accuracy, backward compatibility, and reproducibility.

## Scientific conventions

- Preserve the normalization conventions used throughout the code, including the distinction between gyroradius and inertial-length normalizations.
- Maintain consistent definitions of frequencies, wavevectors, growth rates, and species parameters.
- Cold species are identified by tauS=0. Preserve the distinction between cold-fluid and kinetic contributions.
- Do not change physical assumptions, equations, or numerical approximations without explaining and validating the change.
- Distinguish physical limitations from numerical limitations.

## Fortran development

- Follow the existing Fortran conventions and module structure.
- Prefer targeted modifications over unnecessary refactoring.
- Preserve public interfaces and existing input/output formats unless explicitly authorized.
- Avoid introducing unnecessary dependencies.
- Document nontrivial numerical algorithms and their assumptions.

## Numerical methods

- Preserve established convergence properties and limiting behavior.
- Pay particular attention to small and large wavevector limits, resonances, and transitions between kinetic and cold-plasma descriptions.
- Preserve legacy numerical implementations where they are maintained for compatibility.
- Verify numerical changes against analytical limits, existing reference results, or both.

## Testing and verification

- Identify and run the relevant tests for any code changes.
- Preserve existing regression results unless a justified numerical correction requires updating them.
- Add regression tests for new physical capabilities or corrected numerical behavior.
- Report tests that were not run or did not pass.
- Never change reference outputs merely to obtain passing tests without identifying the physical or numerical reason.

## Documentation

- Documentation must describe actual implemented behavior.
- Preserve equations, normalization conventions, references, and technical qualifications.
- Avoid promotional language, unnecessary repetition, and unsupported claims.
- Use the scientific-editor skill for substantive prose editing.
- Follow additional instructions in docs/AGENTS.md where applicable.

## Working practices

- Inspect relevant existing implementations before proposing changes.
- Make the smallest changes consistent with the task.
- Do not perform unrelated refactoring.
- Review the final diff for unintended changes.
- Summarize substantive modifications, validation performed, and unresolved issues.
