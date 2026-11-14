# Copy a model

Models are modified in place, so assigning a model to a new variable
does not copy it. `sim_copy()` returns an independent copy.

## Usage

``` r
sim_copy(model, name = model$name)
```

## Arguments

- model:

  A
  [`sim_model()`](https://adamrsophie.github.io/simcompR/reference/sim_model.md).

- name:

  Optional new name.

## Value

A new `sim_model`.

## Examples

``` r
m <- sim_model("a")
m2 <- sim_copy(m, "b")
sim_add_block(m2, "constant", "c", value = 1)
length(m$blocks)
#> [1] 0
length(m2$blocks)
#> [1] 1
```
