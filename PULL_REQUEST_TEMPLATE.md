# NA

## What this does

Closes

## Type

- `feat` — new exported function / capability
- `fix` — bug fix
- `test` — tests only
- `docs` — documentation
- `vignette` — worked analysis
- `refactor` / `chore`

## Checklist

Full roxygen2 docs on any new/changed exported function (`@param`,
`@returns`, `@examples`, `@export`)

Each new exported function belongs to one step: build, intersect or
summarise (see SPEC.md, “Package design”); the PR description says which

Ran `devtools::document()` — `man/` and `NAMESPACE` are up to date

Added/updated tests in the mirroring `tests/testthat/test-*.R`; they use
synthetic data (no external files, keys, or network)

`devtools::check()` passes locally with no new errors/warnings/NOTEs

Bumped the dev version in `DESCRIPTION` if I changed an exported
function

If I added a dependency, I ran `renv::snapshot()` and explain the need
below

## Notes for reviewer
