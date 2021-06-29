message(STATUS "generating benchmarks")

foreach(test_category IN LISTS qe_test_categories)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${test_category})
    message(VERBOSE "generating benchmarks in ${test_category} category")
    set(check_pseudo_test_name "benchmark--${test_category}-pseudo")
    download_pseudo(${check_pseudo_test_name} "${test_category}" ${CMAKE_CURRENT_SOURCE_DIR} "benchmark--${test_category}")

    file(GLOB test_input_file_names
         RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/${test_category}"
         "${test_category}/*.in")

    foreach(test_input_file IN LISTS test_input_file_names)
      string(REGEX REPLACE ".in$" "" test_input_file_name_no_extension ${test_input_file})
      set(test_work_dir "${CMAKE_CURRENT_BINARY_DIR}/${test_category}_${test_input_file_name_no_extension}")
      message(VERBOSE "Test input file name without extension ${test_input_file_name_no_extension}")
      message(VERBOSE "Work directory ${test_work_dir}")

      file(MAKE_DIRECTORY "${test_work_dir}")
      file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/${test_category}/${test_input_file}
                       ${test_work_dir}/${test_input_file} SYMBOLIC)

      set(test_name "benchmark--${test_category}_${test_input_file_name_no_extension}")
      set(test_binary "${qe_BINARY_DIR}/bin/${test_category}.x")

      message(VERBOSE "test name : ${test_name}")
      qe_runner(${test_name} ${test_category} ${test_binary} ${test_input_file} "${test_input_file_name_no_extension}.out")

      set_tests_properties(${test_name} PROPERTIES
                           WORKING_DIRECTORY ${test_work_dir}
                           ENVIRONMENT "ESPRESSO_PSEUDO=${ESPRESSO_PSEUDO};OMP_NUM_THREADS=${actual_test_num_threads}"
                           PROCESSORS ${TOT_PROCS} PROCESSOR_AFFINITY TRUE
                           REQUIRED_FILES ${test_input_file}
                           PASS_REGULAR_EXPRESSION "JOB DONE"
                           LABELS "benchmark--${test_category}"
                           DEPENDS ${check_pseudo_test_name})
    endforeach()
  endif()
endforeach()
