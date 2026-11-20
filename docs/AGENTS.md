
# PLUME documentation

PLUME is a scientific plasma dispersion solver. Documentation must accurately describe implemented physical models, numerical methods, and software behavior.

## Editorial constraints

- Preserve the distinction between implemented capabilities, limitations, and planned work.
- Do not alter mathematical definitions, normalization conventions, equations, or physical claims without explicit authorization.
- Preserve command syntax, examples, references, and links.
- Do not modify Fortran source code during prose-only editing tasks.
- Prefer targeted improvements over wholesale rewriting.
- Retain the author's formal, economical scientific style.

## Documentation structure

- README.md is a standalone introduction for GitHub readers.
- docs/index.md is the documentation homepage, providing orientation and navigation.
- Technical manual pages contain the detailed scientific and numerical explanations.
- Avoid unnecessary repetition between these components while retaining sufficient context for each to stand independently.

## Validation

For documentation changes, inspect the Git diff.
When relevant, run scripts/build_docs.sh to check documentation generation and internal links.
Report scientific ambiguities separately rather than resolving them through editorial changes.
