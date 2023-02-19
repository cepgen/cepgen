find_package(Doxygen)
find_package(Sphinx)
if(NOT DOXYGEN_FOUND)
  return()
endif()

if(SPHINX_FOUND)
  set(DOXYGEN_XML YES)
  set(DOXYGEN_HTML NO)
  set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
else()
  set(DOXYGEN_XML NO)
  set(DOXYGEN_HTML YES)
  set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/docs)
endif()
set(CEPGEN_VERSION "${VERSION}")
string(TIMESTAMP CURRENT_YEAR "%Y")
set(DOXYGEN_INPUT_DIR "${PROJECT_SOURCE_DIR}/CepGen ${PROJECT_SOURCE_DIR}/CepGenAddOns ${PROJECT_SOURCE_DIR}/CepGenProcesses")
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/xml/index.xml)
set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
set(SPHINX_IN ${CMAKE_CURRENT_SOURCE_DIR}/conf.in.py)
set(SPHINX_OUT ${CMAKE_CURRENT_SOURCE_DIR}/conf.py)
configure_file(${SPHINX_IN} ${SPHINX_OUT} @ONLY)
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)
# Doxygen won't create this for us
file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})
message(STATUS "Doxygen build started")
message(STATUS ${DOXYGEN_INDEX_FILE})
add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
                  COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
                  MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
                  COMMENT "Generating documentation with Doxygen"
                  VERBATIM)
add_custom_target(Doxygen DEPENDS ${DOXYGEN_INDEX_FILE})
if(SPHINX_FOUND)
  message(STATUS "Sphinx found in ${SPHINX_EXECUTABLE}")
  set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
  set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/output_html)
  set(SPHINX_INDEX_FILE ${SPHINX_BUILD}/index.html)
  add_custom_target(Sphinx
                    COMMAND ${SPHINX_EXECUTABLE} -b html
                    # Tell Breathe where to find the Doxygen output
                    -Dbreathe_projects.CepGen=${DOXYGEN_OUTPUT_DIR}/xml
                    ${SPHINX_SOURCE} ${SPHINX_BUILD}
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/index.rst
                    ${DOXYGEN_INDEX_FILE}
                    #MAIN_DEPENDENCY ${SPHINX_SOURCE}/conf.py
                    COMMENT "Generating documentation with Sphinx")
endif()

