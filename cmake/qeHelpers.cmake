###########################################################
# QE build helper functions
###########################################################

if(NOT TARGET QEGlobalCompileDefinitions)
    add_library(QEGlobalCompileDefinitions INTERFACE)
endif()

function(qe_add_global_compile_definitions DEF)
    if(TARGET QEGlobalCompileDefinitions)
        set_property(TARGET QEGlobalCompileDefinitions APPEND
                     PROPERTY INTERFACE_COMPILE_DEFINITIONS ${DEF} ${ARGN})
    endif()
endfunction(qe_add_global_compile_definitions)

function(qe_get_global_compile_definitions OUTVAR)
    if(TARGET QEGlobalCompileDefinitions)
        get_target_property(defs QEGlobalCompileDefinitions
            INTERFACE_COMPILE_DEFINITIONS)
        set(${OUTVAR} ${defs} PARENT_SCOPE)
    endif()
endfunction(qe_get_global_compile_definitions)

function(qe_get_fortran_cpp_flag OUTVAR)
    if(DEFINED Fortran_PREPROCESSOR_FLAGS)
        set(${OUTVAR} "${Fortran_PREPROCESSOR_FLAGS}" PARENT_SCOPE)
    else()
        # TODO actual flag check
	if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel") 
	  set (${OUTVAR} "-fpp;-allow;nofpp_comments" PARENT_SCOPE)
        else() 
          set(${OUTVAR} "-cpp" PARENT_SCOPE)
        endif() 
    endif()
endfunction(qe_get_fortran_cpp_flag)

function(qe_preprocess_source IN OUT)
    qe_get_global_compile_definitions(global_defs)
    foreach(DEF ${global_defs})
        list(APPEND global_flags "-D${DEF}")
    endforeach()
    get_filename_component(out_dir ${OUT} DIRECTORY)
    if(NOT EXISTS ${out_dir})
        file(MAKE_DIRECTORY ${out_dir})
    endif()
    add_custom_command(
        OUTPUT ${OUT}
        COMMAND ${QE_CPP_FULL_PATH} -P ${global_flags} -E ${IN} > ${OUT}
        MAIN_DEPENDENCY ${IN}
        COMMENT "Preprocessing ${IN}"
        VERBATIM)    
endfunction(qe_preprocess_source)

function(qe_fix_fortran_modules TGT)
    set(targets ${TGT} ${ARGN})
    foreach(tgt IN LISTS targets)
        get_target_property(tgt_type ${tgt} TYPE)
        # All of the following target modifications make
        # sense on non-interfaces only
        if(NOT ${tgt_type} STREQUAL "INTERFACE_LIBRARY")
            get_target_property(tgt_module_dir ${tgt} Fortran_MODULE_DIRECTORY)
            # set module path to tgt_binary_dir/mod
            get_target_property(tgt_binary_dir ${tgt} BINARY_DIR)
            set_target_properties(${tgt}
                PROPERTIES
                    Fortran_MODULE_DIRECTORY ${tgt_binary_dir}/mod/${TGT})
            # make module directory available for clients of TGT 
            target_include_directories(${tgt}
                PUBLIC
                    $<BUILD_INTERFACE:${tgt_binary_dir}/mod/${TGT}>
                INTERFACE
                    $<INSTALL_INTERFACE:${QE_INSTALL_Fortran_MODULES}/qe/${TGT}>)
        endif()
    endforeach()
endfunction(qe_fix_fortran_modules)

function(qe_enable_cuda_fortran SRCS)
    if(QE_ENABLE_CUDA)
        foreach(src IN LISTS SRCS)
            set_source_files_properties(${src} 
                PROPERTIES
                    COMPILE_OPTIONS "${QE_CUDA_COMPILE_OPTIONS}")
        endforeach()
    endif()
endfunction(qe_enable_cuda_fortran)

function(_qe_add_cuda_link_flags TGT)
    if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI" OR CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
        get_target_property(target_type ${TGT} TYPE)
        if(target_type STREQUAL "EXECUTABLE")
            target_link_options(${TGT}
                PRIVATE
                    ${QE_CUDA_COMPILE_OPTIONS})
        endif()
    endif()
endfunction(_qe_add_cuda_link_flags)

function(qe_git_submodule_update PATH)
    # validate submodule_commit_hash_records against git database
    get_filename_component(SUBMODULE_NAME ${PATH} NAME)
    get_filename_component(SUBMODULE_UPPER_DIR ${PATH} DIRECTORY)
    set(commit_hash_file ${qe_SOURCE_DIR}/${SUBMODULE_UPPER_DIR}/submodule_commit_hash_records)
    # a submodule hash consistency check
    if(IS_GIT_PROJECT AND EXISTS ${commit_hash_file})
        # Extract submodule commit hash from git repo database
        execute_process(COMMAND ${GIT_EXECUTABLE} ls-tree HEAD ${PATH}
                        OUTPUT_VARIABLE DATABASE_STRING
                        WORKING_DIRECTORY ${qe_SOURCE_DIR})
        string(REGEX REPLACE " |\t" ";" DATABASE_OUTPUT ${DATABASE_STRING})
        list(GET DATABASE_OUTPUT 2 DATABASE_HASH)

        # Extract submodule commit hash from saved records
        execute_process(COMMAND grep ${SUBMODULE_NAME} ${commit_hash_file}
                        OUTPUT_VARIABLE RECORD_STRING)
        string(REPLACE " " ";" RECORD_OUTPUT ${RECORD_STRING})
        list(GET RECORD_OUTPUT 0 RECORD_HASH)

        if(NOT DATABASE_HASH STREQUAL RECORD_HASH)
          message(FATAL_ERROR "If you are a user, please file a bug report! "
                              "If you are a developer, probably submodule '${SUBMODULE_NAME}' is being touched. "
                              "Inconsistent submodule commit hashes have been detected.\n"
                              "  ${DATABASE_HASH} from repo data base.\n"
                              "  ${RECORD_HASH} from 'submodule_commit_hash_records' file.\n")
        endif()
    endif()

    if(IS_GIT_PROJECT)
        # Old versions of git aren't able to run init+update
        # in one go (via 'git submodule update --init'), we need
        # to call one command for each operation:
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule init -- ${PATH}
                        WORKING_DIRECTORY ${qe_SOURCE_DIR})
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --depth 5 -- ${PATH}
                        WORKING_DIRECTORY ${qe_SOURCE_DIR})
    else()
        if(EXISTS ${commit_hash_file})
            if(EXISTS ${qe_SOURCE_DIR}/${PATH}/.git)
                message(STATUS "Previous clone found at ${qe_SOURCE_DIR}/${PATH}.")
            else()
                # get repo URL
                execute_process(COMMAND ${GIT_EXECUTABLE} config --file .gitmodules --get submodule.${PATH}.URL
                                OUTPUT_VARIABLE SUBMODULE_URL
                                WORKING_DIRECTORY ${qe_SOURCE_DIR}
                                OUTPUT_STRIP_TRAILING_WHITESPACE)

                # Extract submodule commit hash from saved records
                execute_process(COMMAND grep ${SUBMODULE_NAME} ${commit_hash_file}
                                OUTPUT_VARIABLE RECORD_STRING)
                string(REPLACE " " ";" RECORD_OUTPUT ${RECORD_STRING})
                list(GET RECORD_OUTPUT 0 RECORD_HASH)

                message(STATUS "Cloning ${SUBMODULE_URL} into ${qe_SOURCE_DIR}/${PATH}.")

                execute_process(COMMAND ${GIT_EXECUTABLE} init ${PATH}
                                WORKING_DIRECTORY ${qe_SOURCE_DIR})
                execute_process(COMMAND ${GIT_EXECUTABLE} remote add origin ${SUBMODULE_URL}
                                WORKING_DIRECTORY ${qe_SOURCE_DIR}/${PATH})
                execute_process(COMMAND ${GIT_EXECUTABLE} fetch --depth 1 origin ${RECORD_HASH}
                                WORKING_DIRECTORY ${qe_SOURCE_DIR}/${PATH}
                                RESULT_VARIABLE GIT_FETCH_FAILED)
                if(GIT_FETCH_FAILED)
                    file(REMOVE_RECURSE ${qe_SOURCE_DIR}/${PATH}/.git)
                    message(FATAL_ERROR "git fetch failed! Be sure to make ${qe_SOURCE_DIR}/${PATH} completely empty "
                                        "(no hidden files) or removed before a re-try. "
                                        "If this was caused by the lack of Internet access, "
                                        "you can execute 'initialize_external_repos.sh' under QE 'TOPDIR/external' subdirectory "
                                        "on a laptop and transfer the whole 'external' subdirectory.")
                endif()

                execute_process(COMMAND ${GIT_EXECUTABLE} checkout -b recorded_HEAD FETCH_HEAD
                                WORKING_DIRECTORY ${qe_SOURCE_DIR}/${PATH})
            endif()
        else()
          message(FATAL_ERROR "Failed to handle submodule '${SUBMODULE_NAME}'!")
        endif()
    endif()
endfunction(qe_git_submodule_update)

function(qe_add_executable EXE)
    add_executable(${EXE} ${ARGN})
    _qe_add_target(${EXE} ${ARGN})
endfunction(qe_add_executable)

function(qe_add_library LIB)
    add_library(${LIB} ${ARGN})
    set_target_properties(${LIB} PROPERTIES
        SOVERSION ${PROJECT_VERSION_MAJOR}
        VERSION ${PROJECT_VERSION}
    )
    _qe_add_target(${LIB} ${ARGN})
endfunction(qe_add_library)

# Only use this one for Fortran targets
function(_qe_add_target TGT)
    if(TARGET QEGlobalCompileDefinitions)
        target_link_libraries(${TGT} PUBLIC QEGlobalCompileDefinitions)
    endif()
    qe_fix_fortran_modules(${TGT})
    qe_get_fortran_cpp_flag(f_cpp_flag)
    target_compile_options(${TGT} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${f_cpp_flag}>)
    if(QE_ENABLE_CUDA)
        _qe_add_cuda_link_flags(${TGT})
    endif()
endfunction(_qe_add_target)

function(qe_install_targets TGT)
    set(targets ${TGT} ${ARGN})
    install(TARGETS ${targets}
        EXPORT qeTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # Windows needs RUNTIME also for libraries
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/qe/${TGT})
    # Retrieving non-whitelisted properties leads to an hard
    # error, let's skip the following section for interface
    # targets. See here for details:
    # https://gitlab.kitware.com/cmake/cmake/issues/17640
    foreach(tgt IN LISTS targets)
        get_target_property(tgt_type ${tgt} TYPE)
        if(NOT ${tgt_type} STREQUAL "INTERFACE_LIBRARY")
            # If the target generates Fortran modules, make sure
            # to install them as well to a proper location
            get_target_property(tgt_module_dir ${tgt} Fortran_MODULE_DIRECTORY)
            if(tgt_module_dir)
                install(DIRECTORY ${tgt_module_dir}/
                    DESTINATION ${QE_INSTALL_Fortran_MODULES}/qe/${TGT})
            endif()
        endif()        
    endforeach()
endfunction(qe_install_targets)

function(qe_ensure_build_type DEFAULT)
    if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
        message(STATUS "Setting build type to '${DEFAULT}' as none was specified")
        set(CMAKE_BUILD_TYPE "${DEFAULT}"
            CACHE STRING "Choose the type of build." FORCE)
        # Set the possible values of build type for cmake-gui
        set_property(CACHE CMAKE_BUILD_TYPE
            PROPERTY
                STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
    endif()
endfunction(qe_ensure_build_type)

if(TARGET QEGlobalCompileDefinitions)
    qe_install_targets(QEGlobalCompileDefinitions)
endif()
