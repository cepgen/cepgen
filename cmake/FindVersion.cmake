find_package(Git)
set(GIT_HASH "N/A")
set(GIT_DIFF "")
set(GIT_TAG "N/A")
set(GIT_BRANCH "N/A")
if(Git_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} log --pretty=format:'%h' -n 1 OUTPUT_VARIABLE GIT_HASH ERROR_QUIET)
    if(NOT "${GIT_HASH}" STREQUAL "")
        execute_process(COMMAND bash -c "${GIT_EXECUTABLE} diff --quiet --exit-code || echo +" OUTPUT_VARIABLE GIT_DIFF)
        execute_process(COMMAND ${GIT_EXECUTABLE} describe --exact-match --tags OUTPUT_VARIABLE GIT_TAG ERROR_QUIET)
        execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD OUTPUT_VARIABLE GIT_BRANCH)
        string(STRIP "${GIT_HASH}" GIT_HASH)
        string(SUBSTRING "${GIT_HASH}" 1 7 GIT_HASH)
        string(STRIP "${GIT_DIFF}" GIT_DIFF)
        string(STRIP "${GIT_TAG}" GIT_TAG)
        string(STRIP "${GIT_BRANCH}" GIT_BRANCH)
    endif()
endif()

set(VERSION
"#include \"CepGen/Version.h\"

namespace cepgen {
  const std::string version::tag = \"${GIT_TAG}\";
  const std::string version::extended = \"${GIT_HASH}${GIT_DIFF}(${GIT_BRANCH})\";
}")

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/Version.cpp)
    file(READ ${CMAKE_CURRENT_SOURCE_DIR}/Version.cpp VERSION_)
else()
    set(VERSION_ "")
endif()

if(NOT "${VERSION}" STREQUAL "${VERSION_}")
    file(WRITE ${CMAKE_CURRENT_SOURCE_DIR}/Version.cpp "${VERSION}")
endif()
