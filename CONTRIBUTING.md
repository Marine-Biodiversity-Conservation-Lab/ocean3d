# Contributing to ocean3d

Thanks for your interest in contributing! **ocean3d** is an R package for
three-dimensional marine spatial analysis. It
grew out of several lab analyses (Bangladesh fisheries overlap, WOA
environmental extraction, deep-sea depth refuge) and is being generalised into
reusable, tested functions.

This document covers how to get set up, how we branch and coordinate, and what
a contribution needs to look like to be merged.

- [Quick orientation](#quick-orientation)
- [Development setup](#development-setup)
- [Data setup](#data-setup)
- [Branching and coordination](#branching-and-coordination)
- [What a good contribution looks like](#what-a-good-contribution-looks-like)
- [Coding conventions](#coding-conventions)
- [Testing](#testing)
- [Documentation and vignettes](#documentation-and-vignettes)
- [Known rough edges / good first issues](#known-rough-edges--good-first-issues)

---

## Why contribute to **ocean3d**? 

**ocean3d** is an open-source package, built on collaboration across researchers and projects. Currently, **ocean3d** is a very new project undergoing active development, which means that you can shape and directly contribute to the work done. Functionality that you wished existed? You can write an issue with what you would like, reaching contributors with your direct feedback.  Even better, you can develop that functionality and add it to **ocean3d** yourself, broadening the impact of your work. Other people can then come along and use your code to use in their projects, making science overall more efficient and impactful. 

---

## Useful resources for learning how to contribute

### Writing code in R
R Packages (2e)
https://r-pkgs.org/

Advanced R
http://adv-r.had.co.nz/ 

### Using Git and GitHub
Happy Git and GitHub for the useR
https://happygitwithr.com/ 

GitHub Flow
https://docs.github.com/en/get-started/using-github/github-flow

---

## How contributing works for **ocean3d** 

### The loop
```
GitHub issue  →  self-assign  →  branch  →  PR  →  merge
```

We use issues as the unit of work and the record of its current state. 

1. **Find or open an issue.** Every piece of work starts as an issue. 
    - Issues describe anything that needs to be developed, fixed, refined within this codebase. The easiest way to open issues is via this package's [GitHub repository](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/issues). When creating the issue, choose any template that fits (or write your own from scratch). You can write issues without implementing it; it's totally acceptable and encouraged to write issues about what you wish this package did or problems that you have found with it, without intention of coding those actual fixes yourself. 

2. **Claim it.** 
    - Assign yourself and comment before you start coding. You can claim your own issues, or claim another issue that someone else has written but it currently unassigned. This is how we avoid two people extracting the same function from the same prior analysis. If an issue has been assigned, but it's been silent for a month, comment and see if the assigned contributor is still working on it. 
3. **Branch from an up-to-date `main`.**
    - Create a new branch from the most recent commit on `main`. Easiest way is click `Create a branch` under `Development`, within the sidebar of the Issue page in GitHub. We use short-lived branches off `main`, with **one issue per branch**.
    - Branches will be named with the format of `<issue-number>-<short-slug>`. Using the above `Create a branch` procedure should automatically generate branch names based on the issue number and title of the issue. 
4. **Open a Pull Request early** 
    - A draft Pull Request is fine. This way, others can see what areas are being worked on. A Pull Request is a review of the work done on the branch, so that it can be merged with the main branch. This is the main gate for coordinating work done. Ideally, we would have some formal process for PRs, but for now Jay will review these before merging with main to keep work coordinated.
5. **Get Pull Request approval**
    - Another contributor should review and approve of your code before the PR can get merged into the main branch. 
5. **Merge** and delete the branch. Currently, [@JayMatsushiba](https://github.com/JayMatsushiba) will be in charge of this step. 
    - Reference the issue from the PR so it closes automatically. This should happen automatically if the branch is created from the issue. 

`main` must always be installable — never push work-in-progress directly to it. All changes to `main` should be done via a Pull Request. 

### Ownership map — how to avoid collisions

Work is naturally partitioned by file. Coordinate before crossing these lines,
because these files are where merge conflicts actually happen:

| Area | Files | Typical SPEC section |
| --- | --- | --- |
| Volume / voxelisation | `R/volume.R` | Volume calculation |
| Spatial query | `R/intersect.R` | Spatial query |
| Bathymetry | `R/gebco_bathymetry.R` | Data loading and preparation |
| Environmental extraction | `R/extract_to_area.R`, `R/extract_to_points.R`, `R/woa.R` | Environmental extraction, WOA utilities |
| Fisheries | `R/gfw.R` | Data source utilities (GFW) |
| Species data | `R/iucn_utils.R` | Data loading and preparation |

Two people *can* work in the same area — just make the split explicit in the
issue thread.

### Labels

Please apply these when you open an issue:

The templates apply a type label for you; add priority and area yourself.

- Priority: `high`, `medium`, `low` 
- Type (one per branch type): `bug`, `enhancement`, `documentation`, `test`,
  `refactor`, `chore`
- Extra: `good first issue`, `needs-data` (= blocked on a dataset the
  contributor may not have)

### Review and merge

- Every PR needs one approving review from a maintainer (currently
  [@JayMatsushiba](https://github.com/JayMatsushiba)) before merge.
- CI (`R CMD check`) must be green.
- Merge with a message referencing the issue:
  `Add load_eez() for Marine Regions EEZ polygons (#12)`.
- Delete the branch after merging.

### Versioning

Development version is `0.1.0`. Bump the third component in `DESCRIPTION`
(`0.1.1`, …) when a PR adds or changes an exported function. 

---

## What a good contribution looks like

A PR implementing a SPEC function should include **all** of:

1. The function in the appropriate `R/*.R` file, with full roxygen2 docs
   (`@param`, `@returns`, `@examples`, `@export`).
2. Tests in the mirroring `tests/testthat/test-*.R` file, covering the happy
   path, input validation errors, and at least one edge case.
3. `devtools::document()` run, so `man/` and `NAMESPACE` are updated.
4. `devtools::check()` passing locally with no new NOTEs.

If you are generalising code from one of the prior analyses, the goal is **generalisation**, not a copy-paste.
Parameterise the hard-coded species lists, study areas, and file paths, and make the function work for any input meeting the documented contract.

---

## Quick orientation

| Path | What lives there |
| --- | --- |
| `R/` | Package functions, one file per topic |
| `tests/testthat/` | One test file per `R/` file (`test-volume.R` ↔ `R/volume.R`) |
| `tests/testthat/_vcr/` | Recorded HTTP fixtures (cassettes) for API tests |
| `man/` | roxygen2-generated docs — **never edit by hand** |
| `vignettes` | Shorter examples or tutorials; these are easy to follow directions for using the package. Included in package build |
| `vignettes/articles` | Long-form worked analyses; these are the reproductions of the source papers. Not included in package build |
| `SPEC.md` | Why the package works the way it does, decisions taken, and specs for work not yet built. |
| `renv.lock` | Pinned dependency versions |

**Read `SPEC.md` before proposing work.** Its *Architecture and conventions*
section is binding on any new code, and its *Planned work* section specifies the
things that still need building — each entry is written to seed an issue, with an
intended signature and, where one exists, a note on the prior analysis it should
be generalised from (that source material is not in this repo — ask the
maintainer).

To see what needs working on, see the [open
issues](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/issues)
and the package reference (`man/`, `?ocean3d`).

---

## Development setup

### 1. Requirements

- **R 4.5.3** 
- Git, and system GDAL/GEOS/PROJ + NetCDF libraries for `sf` and `terra`
  - Ubuntu/Debian: `sudo apt install libgdal-dev libgeos-dev libproj-dev libudunits2-dev netcdf-bin libnetcdf-dev`
  - macOS: `brew install gdal geos proj udunits netcdf`

### 2. Clone and restore the environment

```bash
git clone https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d.git
cd ocean3d
```

```r
renv::restore()   # installs pinned dependency versions
devtools::load_all()
```

If `renv::status()` reports the project is out of sync, run `renv::restore()`
to match the lockfile. Only run `renv::snapshot()` — and commit the resulting
`renv.lock` change — when your PR **intentionally** adds or upgrades a
dependency; say so in the PR description.

Note that `gfwr` is installed from GitHub (`Remotes:` field in `DESCRIPTION`),
not CRAN.

### 3. API keys

Several functions and vignettes need credentials. Put them in `~/.Renviron`
(`usethis::edit_r_environ()`) — **never** commit them:

```
IUCN_REDLIST_KEY="..."   # https://api.iucnredlist.org/  — fetch_species_assessments()
GFW_TOKEN="..."          # https://globalfishingwatch.org/our-apis/tokens — gfwr::gfw_auth()
```

Code and tests must degrade gracefully without these: use
`testthat::skip_if(identical(Sys.getenv("IUCN_REDLIST_KEY"), ""))` or a
recorded `vcr` cassette rather than assuming a key exists.

### 4. Everyday commands

```r
devtools::load_all()          # load package for interactive work
devtools::document()          # regenerate man/ + NAMESPACE from roxygen comments
devtools::test()              # run the full test suite
testthat::test_file("tests/testthat/test-volume.R")   # one file
devtools::check()             # full R CMD check — run before opening a PR
pkgcheck::pkgcheck()          # pkgcheck for rOpenSci submission
```

---

## Data setup

| Dataset | Used for | Where to get it |
| --- | --- | --- |
| GEBCO 2025 sub-ice topo (`.nc`) | The common bathymetry grid for everything | <https://www.gebco.net/data_and_products/gridded_bathymetry_data/> |
| IUCN Red List `SHARKS_RAYS_CHIMAERAS` range shapefile | Species ranges | <https://www.iucnredlist.org/resources/spatial-data-download> (requires a data request) |
| Marine Regions World EEZ v12 | Study-area clipping | <https://www.marineregions.org/downloads.php> |
| World Ocean Atlas 2023 | Environmental covariates | Downloaded automatically by `woa_download()` into a cache dir; see `woa_cache_dir()` |
| Global Fishing Watch effort | Fisheries effort | Fetched via `gfwr` with `GFW_TOKEN` |
| Bangladesh participatory-mapping fishery footprints | `bangladesh-fisheries-3d-overlap` vignette | Not public — contact the maintainer |

## Coding conventions

- **Style:** tidyverse style (`styler::style_pkg()` if in doubt). Snake_case
  function and argument names.
- **Imports:** `terra` is imported wholesale; `sf`, `stringr`, `ggplot2`,
  `rlang`, `magrittr` are imported selectively via `@importFrom`. Elsewhere
  prefer explicit `pkg::fun()` calls. Add new dependencies to `DESCRIPTION`
  deliberately — this package is already heavy, so justify each one in the PR.
- **Optional dependencies:** guard with `requireNamespace("pkg", quietly =
  TRUE)` and give a clear error telling the user what to install (see
  `fetch_species_assessments()` and `rredlist`).
- **Errors:** validate inputs early and fail with an actionable message naming
  the offending argument.
- **Depth layers:** anything producing a multi-depth `SpatRaster` **must** emit
  `{variable}_depth={value}` layer names.
- **Units:** volumes in km³, areas in km², depths in metres as **positive
  numbers increasing downward** (bathymetry from GEBCO is negative below sea
  level — be explicit in your docs about which convention a given argument
  uses).
- **No side effects on load**, no `library()` calls inside `R/`, no writing to
  the user's filesystem outside `tools::R_user_dir()` (see `woa_cache_dir()`
  for the pattern).

---

## Testing

- `testthat` edition 3. One test file per source file.
- **Tests must not require external data, API keys, or network access.**
  Construct small synthetic objects in the test itself —
  see `tests/testthat/test-volume.R` for the pattern to copy.
- HTTP-dependent tests use `vcr` cassettes stored in `tests/testthat/_vcr/`. If
  you add an API-backed test, record a cassette rather than hitting the live
  service; scrub API keys from the recorded YAML before committing.
- Use `skip_if_not_installed()` for optional-dependency tests and
  `skip_on_cran()` for anything slow.
- Keep the suite fast. If a test needs more than a few seconds, it probably
  needs smaller synthetic inputs.

---

## Documentation and vignettes

- Roxygen2 with markdown enabled. Every exported function needs a runnable
  `@examples` block (wrap in `\dontrun{}` only if it genuinely needs external
  data or a key — prefer building a tiny synthetic example instead).
- Regenerate docs with `devtools::document()`; never hand-edit `man/`.
- `README.md` is generated from `README.Rmd` — edit the `.Rmd` and knit with
  `devtools::build_readme()`.
- Vignettes reproduce real analyses end-to-end. Keep chunks `eval = FALSE`
  unless the analysis can run from data the package can fetch itself. Structure
  them as: study question → data loading → analysis → figures → interpretation.

---

## Questions, and how to get help

- For usage, methods, or interpretation questions, open a
  [GitHub Discussion](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/discussions).
- For bugs or feature/function work, open a
  [GitHub issue](https://github.com/Marine-Biodiversity-Conservation-Lab/ocean3d/issues).
- For access to non-public data (the Bangladesh fishery footprints), contact
  the maintainer directly: Jay Matsushiba <hello@jmatsushiba.com>.

## Attribution

Contributors who add substantive code or analysis will be added to `Authors@R`
in `DESCRIPTION` with an appropriate role (`aut` for substantial
contributions, `ctb` otherwise). If your contribution feeds into a manuscript,
authorship will be discussed openly on the relevant issue. Please add your
ORCID when you're added.

By contributing you agree that your contributions are licensed under the
project's [GNU General Public License v3.0](LICENSE.md).
