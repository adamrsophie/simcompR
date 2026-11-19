# Parameter sweeps and Monte Carlo experiments

`sim_sweep()` runs one simulation per row of `parameters`, a data frame
(or list of named lists) of model parameter values. The model is
compiled once; each run uses
[`sim_set_params()`](https://adamrsophie.github.io/simcompR/reference/sim_set_params.md),
which reuses the dependency analysis whenever the structure is
unchanged.

## Usage

``` r
sim_sweep(
  x,
  parameters,
  tspan,
  times = NULL,
  summary = NULL,
  keep = TRUE,
  cores = 1,
  ...
)

sim_monte_carlo(
  x,
  n,
  sampler,
  seed,
  tspan,
  times = NULL,
  summary = NULL,
  keep = TRUE,
  cores = 1,
  ...
)
```

## Arguments

- x:

  A
  [`sim_model()`](https://adamrsophie.github.io/simcompR/reference/sim_model.md)
  or a compiled `sim_plan`.

- parameters:

  Data frame of parameter values, one run per row, or a list of named
  lists.

- tspan, times:

  Passed to
  [`sim_run()`](https://adamrsophie.github.io/simcompR/reference/sim_run.md).

- summary:

  Optional function of a `sim_result` returning a named numeric vector
  with unique, nonempty names. Values are matched to columns by name,
  never by position: a run that returns the names in another order keeps
  their meaning, a run that omits a name gets `NA` there, and a new name
  adds a column (`NA` for the other runs). A run whose summary is not a
  named numeric vector is marked `"failed"` with the reason in
  `$messages`.

- keep:

  Keep the full results (`TRUE`) or only the summaries.

- cores:

  Number of parallel workers.

- ...:

  Further arguments for
  [`sim_run()`](https://adamrsophie.github.io/simcompR/reference/sim_run.md).

- n:

  Number of Monte Carlo runs.

- sampler:

  Function of the run index returning a named list of parameter values.

- seed:

  Integer seed for `sampler`.

## Value

A `sim_sweep` object with elements `parameters`, `results` (when
`keep = TRUE`), `status`, `messages` and `summary`. `parameters` has one
row per run, a `run` column and one column per parameter that any run
sets; a run that does not set a parameter shows the model's value.
Vector-valued parameters are list columns. `summary` is `parameters`
joined with the summary values and the run status.

## Details

`sim_monte_carlo()` draws `n` parameter sets with `sampler(i)` under a
fixed `seed` and runs them like a sweep. The sampler runs with the
Mersenne-Twister, Inversion and Rejection generators whatever the
caller's [`RNGkind()`](https://rdrr.io/r/base/Random.html), so the draws
are reproducible. Afterwards the caller's `.Random.seed` and generator
kinds are restored. One exception cannot be repaired from R: with
`normal.kind = "Box-Muller"`, R keeps the second deviate of each pair
outside `.Random.seed`, and that cached value is lost, so the caller's
next normal draw differs; a warning says so. `white_noise` blocks do not
use R's generator at all.

With `cores > 1` runs use
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html) (not
available on Windows, where runs are sequential). Results do not depend
on `cores`.

## Examples

``` r
m <- sim_model("decay", params = list(a = 1))
sim_add_block(m, "integrator", "x", initial = 1)
sim_add_block(m, "gain", "k", gain = sim_param("a"))
sim_add_block(m, "gain", "neg", gain = -1)
sim_connect(m, "x", "k")
sim_connect(m, "k", "neg")
sim_connect(m, "neg", "x")
s <- sim_sweep(m, data.frame(a = c(0.5, 1, 2)), tspan = c(0, 1),
  summary = function(r) c(x_end = tail(sim_signal(r, "x"), 1)))
s$summary
#>   run   a     x_end  status
#> 1   1 0.5 0.6065307 success
#> 2   2 1.0 0.3678794 success
#> 3   3 2.0 0.1353353 success
mc <- sim_monte_carlo(m, n = 5, sampler = function(i) list(a = runif(1, 0, 2)),
  seed = 42, tspan = c(0, 1), keep = FALSE,
  summary = function(r) c(x_end = tail(sim_signal(r, "x"), 1)))
mc$summary
#>   run         a     x_end  status
#> 1   1 1.8296121 0.1604758 success
#> 2   2 1.8741508 0.1534852 success
#> 3   3 0.5722791 0.5642380 success
#> 4   4 1.6608953 0.1899688 success
#> 5   5 1.2834910 0.2770683 success
```
