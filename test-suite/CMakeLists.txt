set(qe_test_categories "pw" "cp")

set(ESPRESSO_PSEUDO "${qe_SOURCE_DIR}/pseudo" CACHE STRING "Pseudopotential files")
set(NETWORK_PSEUDO "http://www.quantum-espresso.org/wp-content/uploads/upf_files" CACHE STRING "URL to remote pseudopotential folder")
set(TESTCODE_NPROCS 4 CACHE STRING "Number of MPI processes for each test")
set(TESTCODE_NTHREADS 1 CACHE STRING "Number of OpenMP threads for each test")

#set(ENV{ESPRESSO_PSEUDO} ${ESPRESSO_PSEUDO})
#set(ENV{NETWORK_PSEUDO} ${NETWORK_PSEUDO})

#file(MAKE_DIRECTORY "${qe_BINARY_DIR}/test-suite")

foreach(test_category IN LISTS qe_test_categories)
  message(STATUS "generating tests in ${test_category} category")
  #EXECUTE_PROCESS( COMMAND pwd ECHO_OUTPUT_VARIABLE)
  #EXECUTE_PROCESS( COMMAND "./check_pseudo.sh" "${test_category}_"
  #                 WORKING_DIRECTORY "${qe_SOURCE_DIR}/test-suite" ECHO_OUTPUT_VARIABLE)
  set(check_pseudo_test_name "test-suite_pseudo_${test_category}")
  add_test(NAME ${check_pseudo_test_name}
           COMMAND "./check_pseudo.sh" "${test_category}_"
           WORKING_DIRECTORY "${qe_SOURCE_DIR}/test-suite")
  set_tests_properties(${check_pseudo_test_name} PROPERTIES
                       ENVIRONMENT "ESPRESSO_PSEUDO=${ESPRESSO_PSEUDO};NETWORK_PSEUDO=${NETWORK_PSEUDO}"
                       RESOURCE_LOCK shared_pseudo_folder_lock
                       LABELS "integration;${test_category}")

  file(GLOB category_subfolder_names
       LIST_DIRECTORIES TRUE
       RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
       "${test_category}_*")
  #message(${category_subfolder_names} ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
  foreach(subfolder_name IN LISTS category_subfolder_names)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${subfolder_name}")
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                    ${CMAKE_CURRENT_SOURCE_DIR}/${subfolder_name} input_and_reference 
                    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${subfolder_name}")
    file(GLOB test_input_file_names
       RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/${subfolder_name}"
       "${subfolder_name}/*.in")
    #message("${subfolder_name}")
    #message("${test_input_file_names}")
    foreach(test_input_file IN LISTS test_input_file_names)
      if (NOT test_input_file MATCHES "^benchmark")
        string(REGEX REPLACE ".in$" "" test_input_file_name_no_extension ${test_input_file})

        set(test_work_dir "${CMAKE_CURRENT_BINARY_DIR}/${subfolder_name}/${test_input_file_name_no_extension}")
        file(MAKE_DIRECTORY "${test_work_dir}")

        set(test_name "integration--${subfolder_name}--${test_input_file_name_no_extension}")
        set(test_binary "${qe_BINARY_DIR}/bin/${test_category}.x")
        set(test_input_relative "../input_and_reference/${test_input_file}")
        math(EXPR TOT_PROCS "${TESTCODE_NPROCS} * ${TESTCODE_NTHREADS}")
        #message("${subfolder_name} on ${TOT_PROCS} processes_x_threads")

        if (QE_EANBLE_MPI)
          add_test(NAME ${test_name}
                   COMMAND ${test_binary} -inp ${test_input_relative}
                   WORKING_DIRECTORY ${test_work_dir})
        else()
          add_test(NAME ${test_name}
                   COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${TESTCODE_NPROCS} ${MPIEXEC_PREFLAGS} ${test_binary} -inp ${test_input_relative}
                   WORKING_DIRECTORY ${test_work_dir})
        endif()
        set_tests_properties(${test_name} PROPERTIES
                             ENVIRONMENT "ESPRESSO_PSEUDO=${ESPRESSO_PSEUDO};OMP_NUM_THREADS=${TESTCODE_NTHREADS}"
                             PROCESSORS ${TOT_PROCS} PROCESSOR_AFFINITY TRUE
                             DEPENDS ${check_pseudo_test_name})
        set_property(TEST ${test_name} APPEND PROPERTY LABELS "integration;${subfolder_name}")
      endif()
    endforeach()
  endforeach()
endforeach()
