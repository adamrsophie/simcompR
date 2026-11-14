# Benchmarks

Timed scripts kept out of the package build (see `.Rbuildignore`) and out
of CRAN checks. Run them against an installed copy of simcompR, because
`pkgload::load_all()` code is not byte-compiled and is slower:

```sh
Rscript benchmarks/run.R
```

Each script prints its elapsed time and the solver statistics. The cost of
a run is roughly proportional to the number of segments (sample hits,
breakpoints, events and kinks) times the cost of one `deSolve::ode()` call,
plus the number of right-hand side evaluations times the number of blocks.
