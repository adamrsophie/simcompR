# Saving and loading models

``` r

library(simcompR)
```

Models can be written to declarative JSON files and read back. A file
describes block types, names, parameters, connections, events and nested
subsystems. It never contains R code that is executed when the file is
read.

## Round trip

``` r

m <- sim_model("feedback", params = list(k = 2))
sim_add_block(m, "step", "r", time = 1)
sim_add_block(m, "subtract", "e")
sim_add_block(m, "gain", "g", gain = sim_param("k"))
sim_add_block(m, "transfer_function", "plant", num = 1, den = c(1, 0.3))
sim_connect(m, "r", "e", input = 1)
sim_connect(m, "plant", "e", input = 2)
sim_connect(m, "e", "g")
sim_connect(m, "g", "plant")

f <- tempfile(fileext = ".json")
sim_write_json(m, f)
m2 <- sim_read_json(f)

r1 <- sim_run(sim_compile(m), c(0, 5))
r2 <- sim_run(sim_compile(m2), c(0, 5))
identical(r1$signals, r2$signals)
#> [1] TRUE
```

The results are bit-for-bit identical because numbers are stored with 17
significant digits and matrices keep their dimensions.

## The format

``` r

cat(head(readLines(f), 40), sep = "\n")
#> {
#>   "schema": "simcompR.model",
#>   "schema_version": 1,
#>   "name": "feedback",
#>   "params": {
#>     "k": {
#>       "$num": [
#>         "2"
#>       ]
#>     }
#>   },
#>   "blocks": [
#>     {
#>       "name": "r",
#>       "type": "step",
#>       "params": {
#>         "time": {
#>           "$num": [
#>             "1"
#>           ]
#>         }
#>       }
#>     },
#>     {
#>       "name": "e",
#>       "type": "subtract",
#>       "params": {}
#>     },
#>     {
#>       "name": "g",
#>       "type": "gain",
#>       "params": {
#>         "gain": {
#>           "$param": "k"
#>         }
#>       }
#>     },
#>     {
#>       "name": "plant",
#>       "type": "transfer_function",
```

The top-level fields are `schema` (`"simcompR.model"`), `schema_version`
(currently 1), `name`, `params`, `blocks`, `connections` and `events`.
Encoded values:

| JSON form | R value |
|----|----|
| `{"$num": ["0.10000000000000001"]}` | numeric vector, exact |
| `{"$matrix": {"nrow": 2, "ncol": 2, "data": [...]}}` | numeric matrix |
| `{"$chr": [...]}`, `{"$lgl": [...]}`, `{"$list": [...]}` | other vectors and lists |
| `{"$param": "k"}` | `sim_param("k")` |
| `{"$model": {...}}` | nested subsystem |
| `{"$transition": {...}}` | [`sim_transition()`](https://adamrsophie.github.io/simcompR/reference/sim_transition.md) |

Files from a newer schema version are refused with a clear error rather
than read incorrectly.

## What can be saved

Everything that is declarative is portable:

- all built-in blocks except `fcn`, `ode` and `signal_function`;
- `expression` blocks, event resets given as strings, and state machine
  guards and actions given as strings;
- subsystems and their instance parameters;
- user block types registered with `portable = TRUE`.

R functions cannot be represented safely. A model that contains them is
refused, with the offending block named, instead of being saved without
them:

``` r

bad <- sim_model("closure")
sim_add_block(bad, "constant", "u")
sim_add_block(bad, "fcn", "f", fun = function(u, t) u^2)
sim_connect(bad, "u", "f")
sim_write_json(bad, tempfile())
#> Error:
#> ! Block 'f' has type 'fcn', which is defined by R functions and cannot be written to a portable model file.
```

Replace such blocks with `expression` blocks where possible:

``` r

ok <- sim_model("portable")
sim_add_block(ok, "constant", "u", value = 3)
sim_add_block(ok, "expression", "f", expr = "u[1]^2")
sim_connect(ok, "u", "f")
sim_write_json(ok, f2 <- tempfile(fileext = ".json"))
sim_signal(sim_run(sim_compile(sim_read_json(f2)), c(0, 1), times = 1), "f")
#> f 
#> 9
```

## Safety

Expression strings are parsed when the model is read and checked against
a whitelist of arithmetic, comparison, logical and elementary functions.
A file that tries to call anything else is rejected before anything
runs:

``` r

txt <- readLines(f2)
txt <- sub("u\\[1\\]\\^2", "system('echo hi')", txt)
writeLines(txt, f3 <- tempfile(fileext = ".json"))
sim_read_json(f3)
#> Error:
#> ! Block 'f' (expression): `expr` uses functions that are not allowed in portable expressions: system.
```

## Lists

[`sim_model_to_list()`](https://adamrsophie.github.io/simcompR/reference/sim_write_json.md)
and
[`sim_model_from_list()`](https://adamrsophie.github.io/simcompR/reference/sim_write_json.md)
give the same representation as an R list, for storing models in other
formats.
