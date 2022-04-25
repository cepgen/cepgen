macro(add_version_definition core_sources)
    find_package(Git)
    set(GIT_HASH "N/A")
    set(GIT_BRANCH "N/A")
    if(Git_FOUND)
        execute_process(COMMAND ${GIT_EXECUTABLE} log --pretty=format:'%h' -n 1 OUTPUT_VARIABLE GIT_HASH ERROR_QUIET)
        if(NOT "${GIT_HASH}" STREQUAL "")
            execute_process(COMMAND ${GIT_EXECUTABLE} describe --exact-match --tags OUTPUT_VARIABLE GIT_TAG ERROR_QUIET)
            execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD OUTPUT_VARIABLE GIT_BRANCH)
            if(NOT "${GIT_TAG}" STREQUAL "" AND NOT "${VERSION}" STREQUAL "${GIT_TAG}")
                set(VERSION ${GIT_TAG} PARENT_SCOPE)
            endif()
            string(STRIP "${GIT_HASH}" GIT_HASH)
            string(SUBSTRING "${GIT_HASH}" 1 7 GIT_HASH)
            string(STRIP "${VERSION}" VERSION)
            string(STRIP "${GIT_BRANCH}" GIT_BRANCH)
        endif()
    endif()

    string(TIMESTAMP YEAR "%Y")
    set(VERSION_FILE
"#include \"CepGen/Version.h\"

namespace cepgen {
  const std::string version::tag = \"${VERSION}\";
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
