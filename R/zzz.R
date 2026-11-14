.onLoad <- function(libname, pkgname) {
  register_source_blocks()
  register_math_blocks()
  register_dynamic_blocks()
  register_logic_blocks()
  register_subsystem_blocks()
  register_state_machine_block()
}
