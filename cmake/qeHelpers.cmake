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
        set(${OUTVAR} "-cpp" PARENT_SCOPE)
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
        COMMAND cpp -P ${global_flags} -E ${IN} > ${OUT}
        MAIN_DEPENDENCY ${IN}
        COMMENT "Preprocessing ${IN}"
        VERBATIM)    
endfunction(qe_preprocess_source)

function(qe_fix_fortran_modules LIB)
    set(targets ${LIB} ${ARGN})
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
                    Fortran_MODULE_DIRECTORY ${tgt_binary_dir}/mod/${LIB})
            # make module directory available for clients of LIB 
            target_include_directories(${tgt}
                PUBLIC
                    $<BUILD_INTERFACE:${tgt_binary_dir}/mod/${LIB}>
                INTERFACE
                    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/qe/${LIB}>)
        endif()
    endforeach()
endfunction(qe_fix_fortran_modules)

function(qe_git_submodule_update PATH)
    find_package(Git)
    # Old versions of git aren't able to run init+update
    # in one go (via 'git submodule update --init'), we need
    # to call one command for each operation:
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule init -- ${PATH}
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update -- ${PATH}
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
endfunction(qe_git_submodule_update)

function(qe_add_executable EXE)
    add_executable(${EXE} ${ARGN})
    _qe_add_target(${EXE} ${ARGN})
endfunction(qe_add_executable)

function(qe_add_library LIB)
    add_library(${LIB} ${ARGN})
    _qe_add_target(${LIB} ${ARGN})
endfunction(qe_add_library)

function(_qe_add_target TGT)
    if(TARGET QEGlobalCompileDefinitions)
        target_link_libraries(${TGT} PUBLIC QEGlobalCompileDefinitions)
    endif()
    qe_fix_fortran_modules(${TGT})
    qe_get_fortran_cpp_flag(f_cpp_flag)
    target_compile_options(${TGT} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${f_cpp_flag}>)
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
                    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/qe/${TGT})
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
