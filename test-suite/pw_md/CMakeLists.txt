# prepare files before testing.
set(${test_name_prefix}--preparation-CREATED TRUE PARENT_SCOPE)
add_test(NAME ${test_name_prefix}--preparation
COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/md_restart_verlet_original.md ${CMAKE_CURRENT_BINARY_DIR}/md_restart_verlet/md_restart_verlet.md)
