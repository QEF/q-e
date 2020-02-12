###########################################################
# QE build framework
# Please use the following functions in place of the
# corresponding CMake builtin
###########################################################

mark_as_advanced(QE_MANDATORY_TARGETS)

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

function(qe_add_executable EXE)
    add_executable(${EXE} ${ARGN})
    target_link_libraries(${EXE} PUBLIC ${QE_MANDATORY_TARGETS})
    qe_fix_fortran_modules(${EXE})
endfunction(qe_add_executable)

function(qe_add_library LIB)
    add_library(${LIB} ${ARGN})
    target_link_libraries(${LIB} PUBLIC ${QE_MANDATORY_TARGETS})
    qe_fix_fortran_modules(${LIB})
endfunction(qe_add_library)

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

function(qe_add_fortran_preprocessor)
    foreach(src IN LISTS ARGN)
        if(CMAKE_Fortran_COMPILER_ID STREQUAL "AppleClang" OR
           CMAKE_Fortran_COMPILER_ID STREQUAL "ARMCC" OR
           CMAKE_Fortran_COMPILER_ID STREQUAL "ARMClang" OR
           CMAKE_Fortran_COMPILER_ID STREQUAL "Clang" OR
           CMAKE_Fortran_COMPILER_ID STREQUAL "Cray" OR
           CMAKE_Fortran_COMPILER_ID STREQUAL "Flang" OR
           CMAKE_Fortran_COMPILER_ID STREQUAL "G95" OR
           CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" OR
           CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR
           CMAKE_Fortran_COMPILER_ID STREQUAL "MSVC" OR
           CMAKE_Fortran_COMPILER_ID STREQUAL "XL" OR
           CMAKE_Fortran_COMPILER_ID STREQUAL "XLClang"
        )
            set_source_files_properties(${src} PROPERTIES COMPILE_FLAGS "-cpp")
        elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
            set_source_files_properties(${src} PROPERTIES COMPILE_FLAGS "-Mpreprocess")
        endif()
    endforeach()
endfunction(qe_add_fortran_preprocessor)