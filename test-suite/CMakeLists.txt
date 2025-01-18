set(qe_test_categories "pw" "cp") # "ph" "epw" "tddfpt" "hp" tests are not ready
if(QE_ENABLE_OSCDFT)
  list(APPEND qe_test_categories "oscdft")
endif()

set(ESPRESSO_PSEUDO "${qe_SOURCE_DIR}/pseudo" CACHE STRING "Pseudopotential files")
set(NETWORK_PSEUDO "https://pseudopotentials.quantum-espresso.org/upf_files" CACHE STRING "URL to remote pseudopotential folder")
set(TESTCODE_NPROCS 4 CACHE STRING "Number of MPI processes for each test")
set(TESTCODE_NTHREADS 3 CACHE STRING "Number of OpenMP threads for each test")

message("\n"
        "Only pw and cp results from ctest are reliable, we are working on making the rest tests work reliably with ctest. "
        "To run non-pw/cp tests, make a softlink of the bin directory to the root of QE source tree and run tests in the test-suite directory under that root.\n")
find_program(BASH_PROGRAM bash)
if(NOT BASH_PROGRAM)
  message(WARNING "bash not found! Tests under ${CMAKE_CURRENT_SOURCE_DIR} not added.")
  return()
endif()

if(QE_ENABLE_MPI)
  set(actual_test_num_procs ${TESTCODE_NPROCS})
else()
  set(actual_test_num_procs 1)
endif()

if(QE_ENABLE_OPENMP)
  set(actual_test_num_threads ${TESTCODE_NTHREADS})
else()
  set(actual_test_num_threads 1)
endif()

math(EXPR TOT_PROCS "${actual_test_num_procs} * ${actual_test_num_threads}")
message(VERBOSE "Each test runs on ${TOT_PROCS} processes_x_threads")

function(qe_runner test_name label test_binary test_input test_output)
  if(QE_ENABLE_MPI)
    add_test(NAME ${test_name}
             COMMAND ${BASH_PROGRAM} ${qe_SOURCE_DIR}/test-suite/ctest_runner.sh ${label} ${test_output} ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${TESTCODE_NPROCS} ${MPIEXEC_PREFLAGS} ${test_binary} -inp ${test_input})
  else()
    add_test(NAME ${test_name}
             COMMAND ${BASH_PROGRAM} ${qe_SOURCE_DIR}/test-suite/ctest_runner.sh ${label} ${test_output} ${test_binary} -inp ${test_input})
  endif()
endfunction()

function(download_pseudo test_name input_file_dir_prefix work_dir labels)
  add_test(NAME ${test_name}
           COMMAND "${qe_SOURCE_DIR}/test-suite/check_pseudo.sh" "${input_file_dir_prefix}"
           WORKING_DIRECTORY "${work_dir}")
  set_tests_properties(${test_name} PROPERTIES
                       ENVIRONMENT "ESPRESSO_PSEUDO=${ESPRESSO_PSEUDO};NETWORK_PSEUDO=${NETWORK_PSEUDO}"
                       RESOURCE_LOCK shared_pseudo_folder_lock
                       LABELS "${labels};pseudo")
endfunction()

# skip checks for disabled features
if(NOT QE_ENABLE_LIBXC)
  list(APPEND CHECK_SKIP_MESSAGE "libxc needed for this functional")
endif()
if (QE_ENABLE_CUDA) 
  list(APPEND CHECK_SKIP_MESSAGE "not ported to GPU")
endif()
# prepare jobconfig and userconfig for testcode.py
file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/jobconfig
                 ${CMAKE_CURRENT_BINARY_DIR}/jobconfig
                 SYMBOLIC COPY_ON_ERROR)
string(REPLACE "/" "\\/" qe_SOURCE_DIR_FOR_SED ${qe_SOURCE_DIR})
execute_process(COMMAND sed "s/XXXXXX/${qe_SOURCE_DIR_FOR_SED}/" ${CMAKE_CURRENT_SOURCE_DIR}/userconfig.tmp
                OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/userconfig)

foreach(test_category IN LISTS qe_test_categories)
  message(STATUS "generating tests in ${test_category} category")
  set(check_pseudo_test_name "system--${test_category}-pseudo")
  download_pseudo(${check_pseudo_test_name} "${test_category}_" ${CMAKE_CURRENT_SOURCE_DIR} "system--${test_category}")

  file(GLOB category_subfolder_names
       LIST_DIRECTORIES TRUE
       RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
       "${test_category}_*")
  #message(${category_subfolder_names} ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
  foreach(subfolder_name IN LISTS category_subfolder_names)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${subfolder_name}")
    file(GLOB test_input_file_names
       RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/${subfolder_name}"
       "${subfolder_name}/*.in")
    #message("${subfolder_name}")
    #message("${test_input_file_names}")

    # add a correctness check test by leveraging existing testcode.py
    set(test_name_prefix "system--${subfolder_name}")
    set(test_work_dir_prefix "${CMAKE_CURRENT_BINARY_DIR}/${subfolder_name}")
    add_test(NAME ${test_name_prefix}-correctness
                 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/testcode/bin/testcode.py -v -c ${subfolder_name} compare)
    set_tests_properties(${test_name_prefix}-correctness PROPERTIES
                         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                         LABELS "system--${test_category};${subfolder_name}")

    # an adhoc fix for matdyn input missing error from testcode.py
    file(GLOB matdyn_input_file_names
       RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/${subfolder_name}"
       "${subfolder_name}/matdyn.in.*")
    foreach(matdyn_input IN LISTS matdyn_input_file_names)
        file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/${subfolder_name}/${matdyn_input}
                         ${CMAKE_CURRENT_BINARY_DIR}/${subfolder_name}/${matdyn_input} SYMBOLIC)
    endforeach()

    # handle additional test customization
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${subfolder_name}/CMakeLists.txt)
      add_subdirectory(${subfolder_name})
    endif()

    foreach(test_input_file IN LISTS test_input_file_names)
      # skip reference and test run output files which also has filename extension .in
      if (NOT test_input_file MATCHES "^benchmark" AND NOT test_input_file MATCHES "^test")
        file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/${subfolder_name}/${test_input_file}
                         ${CMAKE_CURRENT_BINARY_DIR}/${subfolder_name}/${test_input_file}
                         SYMBOLIC COPY_ON_ERROR)
        set(TESTSUITE_ARGS 0)
        if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${subfolder_name}/benchmark.out.git.inp=${test_input_file})
          file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/${subfolder_name}/benchmark.out.git.inp=${test_input_file}
                           ${CMAKE_CURRENT_BINARY_DIR}/${subfolder_name}/benchmark.out.git.inp=${test_input_file}
                           SYMBOLIC COPY_ON_ERROR)
        else()
          while (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${subfolder_name}/benchmark.out.git.inp=${test_input_file}.args=${TESTSUITE_ARGS} AND TESTSUITE_ARGS LESS 10)
            math(EXPR TESTSUITE_ARGS "${TESTSUITE_ARGS} + 1")
          endwhile()
          if (TESTSUITE_ARGS EQUAL 10)
            #message("Reference output not found for test ${subfolder_name}/${test_input_file}")
            continue()
          endif()
          file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/${subfolder_name}/benchmark.out.git.inp=${test_input_file}.args=${TESTSUITE_ARGS}
                           ${CMAKE_CURRENT_BINARY_DIR}/${subfolder_name}/benchmark.out.git.inp=${test_input_file}.args=${TESTSUITE_ARGS}
                           SYMBOLIC COPY_ON_ERROR)
        endif()
        string(REGEX REPLACE ".in$" "" test_input_file_name_no_extension ${test_input_file})
        if (TESTSUITE_ARGS EQUAL 0)
          set(test_output_file "../test.out.010121.inp=${test_input_file}")
        else()
          set(test_output_file "../test.out.010121.inp=${test_input_file}.args=${TESTSUITE_ARGS}")
        endif()

        # handle test series dependency.
        if(test_input_file_name_no_extension MATCHES "-[1-9]$")
          set(in_test_series TRUE)
          string(REGEX REPLACE "-[1-9]$" "" test_input_file_name_no_extension_no_series ${test_input_file_name_no_extension})

          string(LENGTH ${test_input_file_name_no_extension} name_length)
          math(EXPR last_char_id "${name_length} - 1")
          string(SUBSTRING ${test_input_file_name_no_extension} ${last_char_id} ${last_char_id} last_char)
          #message("----- last char : ${last_char}")
          math(EXPR last_series "${last_char} - 1")

          set(last_test_name "${test_name_prefix}--${test_input_file_name_no_extension_no_series}-${last_series}")
          set(base_test_name "${test_name_prefix}--${test_input_file_name_no_extension_no_series}")
          set(test_work_dir "${test_input_file_name_no_extension_no_series}")
        else()
          set(in_test_series FALSE)
          set(last_test_name "")
          set(base_test_name "")
          set(test_work_dir "${test_input_file_name_no_extension}")
        endif()
        file(MAKE_DIRECTORY "${test_work_dir_prefix}/${test_work_dir}")

        set(test_name "${test_name_prefix}--${test_input_file_name_no_extension}")
        # prefer $<TARGET_FILE:qe_${test_category}_exe>) but currently target names don't match test_category
        set(test_binary "${qe_BINARY_DIR}/bin/${test_category}.x")
        set(test_input_relative "../${test_input_file}")
        #message("test name : ${test_name}")
        #message("last test name : ${last_test_name}")
        #message("base test name : ${base_test_name}")
        #message("test work dir : ${test_work_dir_prefix}/${test_work_dir}")

        qe_runner(${test_name} ${test_category} ${test_binary} ${test_input_relative} ${test_output_file})

        set(expected_exit_msg_file "${CMAKE_CURRENT_SOURCE_DIR}/${subfolder_name}/${test_input_file_name_no_extension}.expected_exit_msg")
        if(EXISTS "${expected_exit_msg_file}")
          file(STRINGS "${expected_exit_msg_file}" expected_exit_msg)
        endif()

        set_tests_properties(${test_name} PROPERTIES
                             WORKING_DIRECTORY ${test_work_dir_prefix}/${test_work_dir}
                             ENVIRONMENT "ESPRESSO_PSEUDO=${ESPRESSO_PSEUDO};OMP_NUM_THREADS=${actual_test_num_threads}"
                             PROCESSORS ${TOT_PROCS} PROCESSOR_AFFINITY TRUE
                             REQUIRED_FILES ${test_input_relative}
                             RESOURCE_LOCK shared_workdir_${test_work_dir}
                             DEPENDS ${check_pseudo_test_name})

        # Each test occupies one GPU regardless of the number of MPI ranks.
        set_tests_properties(${test_name} PROPERTIES RESOURCE_GROUPS "nvidia_gpus:1")

        # capture non-zero exit code
        if(DEFINED expected_exit_msg)
          set_tests_properties(${test_name} PROPERTIES PASS_REGULAR_EXPRESSION "${expected_exit_msg}")
          unset(expected_exit_msg)
        else()
          set_tests_properties(${test_name} PROPERTIES PASS_REGULAR_EXPRESSION "JOB DONE")
        endif()
        set_property(TEST ${test_name} APPEND PROPERTY PASS_REGULAR_EXPRESSION "${CHECK_SKIP_MESSAGE}")

        set_property(TEST ${test_name} APPEND PROPERTY LABELS "system--${test_category};${subfolder_name}")

        # add dependency between workflow steps
        if(in_test_series)
          set_property(TEST ${test_name} APPEND PROPERTY DEPENDS "${base_test_name};${last_test_name}")
        endif()

        # add dependency on test preparation
        if(${test_name_prefix}--preparation-CREATED)
          set_property(TEST ${test_name} APPEND PROPERTY DEPENDS ${test_name_prefix}--preparation)
        endif()

        # make correctness check depend on the test run
        set_property(TEST ${test_name_prefix}-correctness APPEND PROPERTY DEPENDS "${test_name}")
      endif()
    endforeach()
  endforeach()
endforeach()

if(QE_ENABLE_BENCHMARK)
  add_subdirectory(benchmarks)
endif()
