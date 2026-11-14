## Submission

This is a new submission.

## Test environments

* GitHub Actions: ubuntu-latest (R release, R devel and R oldrel-1),
  macOS-latest (R release), windows-latest (R release and R devel), with
  `R CMD check --as-cran`.
* GitHub Actions: ubuntu-latest (R release and R devel) with the CRAN
  incoming checks enabled and the PDF manual built.
* win-builder: R-devel and R-release.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Adam R Sophie <adamrsophie@gmail.com>'
  New submission

## Notes for the reviewers

* 'Simulink', 'Stateflow', 'bdsim', 'deSolve', 'Python' and 'MATLAB' are
  names of software and are quoted in the Description.
* Examples run in well under five seconds each. The vignettes run
  simulations with deSolve and take a few seconds in total.
* The package does not write outside `tempdir()`, does not change the
  user's options, graphics parameters or random number state (the seeded
  noise blocks restore `.Random.seed` and `RNGkind()` on exit), and uses at
  most two cores in examples, tests and vignettes.
