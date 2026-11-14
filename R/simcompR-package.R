#' simcompR: block-diagram modeling and dynamic system simulation
#'
#' `simcompR` builds block diagrams in code, compiles them into a
#' dependency-aware execution plan and simulates continuous, discrete,
#' multirate and hybrid dynamics.
#'
#' The workflow has three stages:
#'
#' 1. Build a model with [sim_model()], [sim_add_block()] and [sim_connect()].
#' 2. Compile it with [sim_compile()]. Compilation validates ports and
#'    signal sizes, detects algebraic loops and fixes the execution order.
#' 3. Run the plan with [sim_run()] and inspect the [sim_result][sim_run()]
#'    with [sim_signal()], [sim_states()], [sim_events()] and `plot()`.
#'
#' See `vignette("simcompR")` for a first example and
#' [sim_block_types()] for the block library.
#'
#' @keywords internal
"_PACKAGE"
