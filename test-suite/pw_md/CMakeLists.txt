# prepare files before testing.
set(${test_name_prefix}--preparation-CREATED TRUE PARENT_SCOPE)

add_test(NAME ${test_name_prefix}--md_restart_verlet-cleanup
COMMAND ${CMAKE_COMMAND} -E rm -rf ${CMAKE_CURRENT_BINARY_DIR}/md_restart_verlet/md_restart_verlet.save)

add_test(NAME ${test_name_prefix}--md_restart_verlet-init
COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/md_restart_verlet_original.md ${CMAKE_CURRENT_BINARY_DIR}/md_restart_verlet/md_restart_verlet.md)

add_test(NAME ${test_name_prefix}--preparation
COMMAND ${CMAKE_COMMAND} -E true)

set_property(TEST ${test_name_prefix}--preparation APPEND PROPERTY DEPENDS "${test_name_prefix}--restart_verlet-cleanup;${test_name_prefix}--md_restart_verlet-init")
