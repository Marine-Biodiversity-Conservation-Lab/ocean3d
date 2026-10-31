# AGENTS.md

Notes for AI coding agents working on ocean3d. Read CONTRIBUTING.md for
the general contribution process and SPEC.md for the package design;
this file covers conventions agents tend to get wrong. When the workflow
changes, update this file in the same PR.

## Data models

ocean3d currently has two data models for gridded data (see SPEC.md,
“Data model”): netCDF files read with ncdf4
([`extract_to_point()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_point.md),
[`extract_to_buffer()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_buffer.md),
[`extract_to_track()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/extract_to_track.md))
and `SpatVoxel` / `SpatEnvelope` objects built on terra (everything
else). \#99 plans to unify them.

- Check which model a function takes before passing data to it or adding
  to it, and say which one in its `@param` docs.
- Don’t add a third representation or a new conversion between the two
  without checking SPEC.md and \#99 first. New work should fit the
  planned model, not widen the split.
- Each analysis function belongs to one step: build, intersect or
  summarise (SPEC.md, “Package design”). Add new features to one step;
  don’t add another function that does all three in one call.

## Reuse before you write

Build on the functions and data representations the package already has
before adding new ones.

- Before writing a helper, search `R/` for one that already does the job
  (reading bounds, detecting variables, decoding time, summarising a
  depth profile). Use it, or extend it with an argument, instead of
  writing a second version for your feature.
- If an existing helper is close but not quite right, change it in one
  place and update its callers, rather than copying it under a new name.
  \#83 did this with a second bounds reader for buffers, which then had
  to be merged back.
- Use the existing classes (`SpatVoxel`, `SpatEnvelope`, `SpatVolume`)
  and their constructors and converters
  ([`as_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_voxel.md),
  [`as_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/as_envelope.md),
  [`vect_to_envelope()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/vect_to_envelope.md),
  [`envelope_to_voxel()`](https://marine-biodiversity-conservation-lab.github.io/ocean3d/reference/envelope_to_voxel.md))
  instead of new structures for the same data.
- Prefer existing dependencies (terra, sf, ncdf4) to hand-written
  geometry or file reading. Check the package does what you need before
  writing your own.
- A new exported function, class or data representation needs a reason
  the existing ones can’t serve. Say what you checked in the PR
  description.

## Internal helpers

Don’t create an internal function unless it earns its place (#97).

Keep a helper when:

- it is called from several places, or
- it holds real logic that is reused or likely to be (e.g. time
  decoding, great-circle distance), or
- it is a clearly named unit with its own tests covering behaviour the
  exported function can’t easily reach.

Inline it when:

- it has one call site and is short,
- it only forwards its arguments to another function or wraps one base R
  call,
- it is an if/else chain that picks a value or an error message for one
  caller,
- it is one layer in a chain of single-use functions.

Large single-use steps (e.g. `.extract_one_file()`) can stay if they
make the exported function read as a sequence of named steps.

When you inline a helper that tests call directly, move those checks
onto the exported function rather than deleting them.

## Conditions (errors, warnings, messages)

- Signal classed conditions with
  [`rlang::abort()`](https://rlang.r-lib.org/reference/abort.html),
  [`rlang::warn()`](https://rlang.r-lib.org/reference/abort.html) and
  [`rlang::inform()`](https://rlang.r-lib.org/reference/abort.html) and
  a `class =` argument. Don’t build conditions by hand with
  `structure(list(message = ..., call = NULL), class = c(...))`.
- Build messages with
  [`glue::glue()`](https://glue.tidyverse.org/reference/glue.html) so
  values appear in place (`"{n} observations have ..."`) instead of long
  [`paste0()`](https://rdrr.io/r/base/paste.html) chains. Work out
  plural wording in a variable first.
- Condition classes are internal: they exist so tests can catch or
  silence a specific condition without matching its wording. Don’t
  document them in help pages, NEWS or README.
- Name classes `ocean3d_<error|warning|message>_<what>`
  (e.g. `ocean3d_warning_timezone`). Some existing classes predate this
  convention; \#97 renames them.
- Attach data to a condition only if package code or a test reads it,
  and give each class a single producer with the same fields every time.
- A helper that signals a condition is worth keeping when the condition
  is raised from several places. With one call site, signal it inline.

## Tests

- Test behaviour through the exported functions, on synthetic data (no
  external files, keys or network).
- Check conditions with `expect_error()` / `expect_warning()` using
  `class =`, plus `expect_match()` on the parts of the message that
  matter.
- Use `skip_if_not_installed()` only for packages in `Suggests` (e.g.
  lwgeom). Packages in `Imports`, such as ncdf4, are always installed,
  so a skip for them never triggers and only adds noise.

## Docs

- roxygen2 docs on every changed exported function. Run
  `devtools::document()`; never edit `man/` or `NAMESPACE` by hand.
- Follow the CONTRIBUTING.md “Coding conventions”, “Testing” and
  “Documentation and vignettes” sections; they aren’t repeated here.

## Workflow

- One agent per issue, on its own branch in its own git worktree,
  branched from up-to-date `origin/main`. Name branches `<N>-<slug>` or
  `<type>/<N>-<slug>` (e.g. `74-coordinate-axis-handling`,
  `docs/59-getting-started-vignette`).
- Put `Closes #N` in the PR body so the issue closes on merge. `(#N)` in
  the title is not enough.
- Keep PRs reviewable in one sitting: about 400 changed lines in `R/`
  and 1000 in total. Generated files (`man/`, `NAMESPACE`, `renv.lock`,
  rendered articles) don’t count. Split bigger work into PRs that each
  branch from `main`, do one thing and pass `devtools::check()`. Open
  the next one after the previous one merges (e.g. shared helpers, then
  the function, then articles). \#83 (a 2100-line `extract_to_buffer.R`)
  is the example to avoid.
- Don’t stack PRs (a PR whose base is another open PR’s branch). Every
  change to a lower PR has to be merged into and re-reviewed in each PR
  above it, which multiplies review work. Branch from `main`, and wait
  for the PR your work depends on to merge.
- Agents can push their branch and open a PR for their issue. A person
  reviews every PR, and only the developer approves and squash-merges:
  don’t merge or approve. Don’t open issues or post comments elsewhere
  on GitHub unless the developer asks.
- Each PR gets an agent code review and a prose review of its docs and
  text. Reply to each review comment with what changed and in which
  commit (e.g. “Fixed in e860cb3: …”).
- Commits and PR descriptions written with an agent include AI
  attribution, e.g. a `Co-Authored-By:` trailer naming the model.
- Before handing back: `devtools::document()`, `devtools::test()`,
  `devtools::check()`, and `lintr::lint()` on changed files with no new
  lints.
