install(DIRECTORY Guib PWgui QE-modes
    DESTINATION
        share/GUI
    USE_SOURCE_PERMISSIONS)

install(FILES README
    DESTINATION
        share/GUI)

install(CODE
    "execute_process(
        COMMAND make init
        WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/share/GUI/PWgui
        OUTPUT_QUIET
        ERROR_QUIET)")
