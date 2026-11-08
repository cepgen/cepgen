macro(cepgen_find_version cepgen_version git_branch git_hash)
    find_package(Git)
    set(${cepgen_version} "N/A")
    set(${git_hash} "N/A")
    set(${git_branch} "N/A")
    if(NOT Git_FOUND)
      return()
    endif()
    execute_process(COMMAND ${GIT_EXECUTABLE} log --pretty=format:'%h' -n 1 OUTPUT_VARIABLE hash ERROR_QUIET)
    string(STRIP "${hash}" hash)
    if(NOT "${hash}" STREQUAL "")
        execute_process(COMMAND ${GIT_EXECUTABLE} describe --first-parent --tags --no-abbrev OUTPUT_VARIABLE git_tag ERROR_QUIET)
        execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD OUTPUT_VARIABLE branch ERROR_QUIET)
        string(STRIP "${git_tag}" git_tag)
        if(NOT "${git_tag}" STREQUAL "")
            set(${cepgen_version} "${git_tag}")
        endif()
        string(SUBSTRING "${hash}" 1 7 ${git_hash})
        string(STRIP "${branch}" ${git_branch})
    endif()
endmacro()

macro(add_version_definition core_sources)
    cepgen_find_version(CEPGEN_VERSION GIT_BRANCH GIT_HASH)
    string(TIMESTAMP YEAR "%Y")
    set(VERSION_FILE
"#include \"CepGen/Version.h\"

namespace cepgen {
  const std::string version::tag = \"${CEPGEN_VERSION}\";
  const std::string version::extended = \"${GIT_HASH}(${GIT_BRANCH})\";
  const std::string version::banner = \"CepGen version \"+version::tag+\" (\"+version::extended+\")\\n\"
    \"Copyright (c) ${YEAR} L. Forthomme.\\n\"
    \"License GPLv3+: GNU GPL version 3 or later <https://gnu.org/licenses/gpl.html>.\\n\"
    \"This is free software: you are free to change and redistribute it.\\n\"
    \"There is NO WARRANTY, to the extent permitted by law.\";
}")
    set(VERSION_FILE_ "")
    set(output "${CMAKE_CURRENT_BINARY_DIR}/Version.cpp")
    if(EXISTS ${output})
        file(READ ${output} VERSION_FILE_)
    endif()
    if(NOT "${VERSION_FILE}" STREQUAL "${VERSION_FILE_}")
        file(WRITE ${output} "${VERSION_FILE}")
    endif()
    list(APPEND core_sources ${output})
endmacro()
