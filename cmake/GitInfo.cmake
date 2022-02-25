set(GITREV_BARE_FILE git-rev.h)
set(GITREV_BARE_TMP git-rev-tmp.h)
set(GITREV_FILE ${qe_BINARY_DIR}/${GITREV_BARE_FILE})
set(GITREV_TMP ${qe_BINARY_DIR}/${GITREV_BARE_TMP})

# Always clean up git-rev-tmp.h file when cmake is run
execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${GITREV_TMP})

# The following custom command picks up changes to the git revision information
# every time the project is rebuilt. Even if the repository is updated (git pull)
# without re-running cmake. It also appends '-dirty' to the commit hash if there are
# unsaved changes to the repository.
#
# To avoid triggering a relink every time, the repository info is saved to
# a temporary file, and the temporary file is copied over the actual file
# only if the contents changed (using 'copy_if_different').
# The temporary file is deleted to force the custom command to run on
# the next build.
#
# Apparently custom commands need to be defined where the output is used.
# If this in the main CMakeLists.txt it does not work.

# Sed flags were once an issue and some HPC have old sed's
set(SED_FLAG "-E")
execute_process(COMMAND "sed" ${SED_FLAG} "s/\"/\\\\\"/g" "<${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt" OUTPUT_QUIET
                ERROR_VARIABLE SED_ERROR)
if(SED_ERROR MATCHES ".*invalid.*")
    set(SED_FLAG "-r")
    execute_process(COMMAND "sed" ${SED_FLAG} "s/\"/\\\\\"/g" "<${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt"
                            OUTPUT_QUIET ERROR_QUIET ERROR_VARIABLE SED_ERROR)
    if(SED_ERROR MATCHES ".*invalid.*")
        message(
            WARNING
                "Your system supports neither the sed -E or -r flag, git revision information will not be included in output"
        )
    else(SED_ERROR MATCHES ".*invalid.*")
        message("   sed supports -r")
    endif(SED_ERROR MATCHES ".*invalid.*")
else(SED_ERROR MATCHES ".*invalid.*")
    message("   sed supports -E")
endif(SED_ERROR MATCHES ".*invalid.*")

add_custom_target(gitrev
  COMMAND ${CMAKE_COMMAND} -E echo_append "#define GIT_BRANCH_RAW \"" > ${GITREV_TMP}
  COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD | sed "s/$/\"/" >> ${GITREV_TMP}
  COMMAND ${CMAKE_COMMAND} -E echo >> ${GITREV_TMP}
  COMMAND ${CMAKE_COMMAND} -E echo_append "#define GIT_HASH_RAW \"" >> ${GITREV_TMP}
  COMMAND ${GIT_EXECUTABLE} describe --always --dirty --abbrev=40 --match="NoTagWithThisName" | sed "s/$/\"/" >> ${GITREV_TMP}
  COMMAND ${CMAKE_COMMAND} -E echo >> ${GITREV_TMP}
  COMMAND ${CMAKE_COMMAND} -E echo_append "#define GIT_COMMIT_LAST_CHANGED_RAW \"" >> ${GITREV_TMP}
  COMMAND ${GIT_EXECUTABLE} log -1 --format=%ad | sed "s/$/\"/" >> ${GITREV_TMP}
  COMMAND ${CMAKE_COMMAND} -E echo >> ${GITREV_TMP}
  COMMAND ${CMAKE_COMMAND} -E echo_append "#define GIT_COMMIT_SUBJECT_RAW \"" >> ${GITREV_TMP}
  COMMAND ${GIT_EXECUTABLE} log -1 --format=%s | cut -c 1-50 | sed ${SED_FLAG} "s/\"/\"\\/\\/'\"'\\/\\/\"/g" | tr -d "\\n" >> ${GITREV_TMP}
  COMMAND ${CMAKE_COMMAND} -E echo_append "\"" >> ${GITREV_TMP}
  COMMAND ${CMAKE_COMMAND} -E echo >> ${GITREV_TMP}
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${GITREV_TMP} ${GITREV_FILE}
  COMMAND ${CMAKE_COMMAND} -E remove ${GITREV_TMP}
  WORKING_DIRECTORY ${qe_SOURCE_DIR}
  VERBATIM)

# Print some configure-time git info (useful for understanding what commits
# are in particular build for the nightly CDash reports)

execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    OUTPUT_VARIABLE GIT_CONFIG_BRANCH
    WORKING_DIRECTORY ${qe_SOURCE_DIR} OUTPUT_STRIP_TRAILING_WHITESPACE)
message("   Git branch: ${GIT_CONFIG_BRANCH}")

# Breaking down the arguments to 'git describe'
#  --abbrev=40     Size of hash to print.  This should print the entire hash.
#  --dirty         Append hash with '-dirty' if there are uncommitted changes.
# The behavior of describe looks for a tag in the parents first, and then falls
# back to the commit hash (if --always is specified)
#  --always        Show the commit hash as fallback
#  --match="NoTagWithThisName"
#     If a tag is found, the output looks like:
#       second_annotated_tag-29-g1fd38cccc0fd2f683ec223ca0783bb671bfedd4e
#     In order to always get just the commit hash, specify a tag pattern
#     that should never match.
execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --always --dirty --abbrev=40 --match="NoTagWithThisName"
    OUTPUT_VARIABLE GIT_CONFIG_COMMIT_HASH
    WORKING_DIRECTORY ${qe_SOURCE_DIR} OUTPUT_STRIP_TRAILING_WHITESPACE)
message("   Git commit hash: ${GIT_CONFIG_COMMIT_HASH}")
