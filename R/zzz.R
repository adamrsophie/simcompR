.onLoad <- function(libname, pkgname) {
  register_source_blocks()
  register_math_blocks()
  register_dynamic_blocks()
  register_logic_blocks()
  register_subsystem_blocks()
  register_state_machine_block()
  # string expressions of these types also see model and subsystem parameters
  for (type in c("expression", "state_machine")) .registry[[type]]$scoped <- TRUE
}
