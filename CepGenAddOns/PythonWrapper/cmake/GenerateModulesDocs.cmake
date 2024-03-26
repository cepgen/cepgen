macro(cepgen_generate_python_descriptions)
    set(options)
    set(one_val NAME)
    set(multi_vals CATEGORY)
    cmake_parse_arguments(ARG "${options}" "${one_val}" "${multi_vals}" ${ARGN})
    exec_program(${MODULE_DOC_EXE}
        ARGS --documentation-generator "\"text<modulesOnly,camelCaseModulesNames\""
             --categories ${ARG_CATEGORY}
             --quiet
        OUTPUT_VARIABLE modules)
    message(STATUS "Python... generating modules description for ${ARG_CATEGORY} in ${ARG_NAME}")
    file(MAKE_DIRECTORY ${MODULE_DOC_OUTPUT}/${ARG_NAME})
    foreach(_mod ${modules})
        exec_program(${MODULE_DOC_EXE}
            ARGS --documentation-generator python
                 --categories ${ARG_CATEGORY}
                 --modules ${_mod}
                 --output "${MODULE_DOC_OUTPUT}/${ARG_NAME}/${_mod}_cfi.py"
                 --quiet
            OUTPUT_VARIABLE dummy)
        message(DEBUG "      ... generated Python configuration file for ${ARG_NAME}/${_mod}\n"
                      "      ... debugging:\n${dummy}")
    endforeach()
endmacro()

file(MAKE_DIRECTORY ${MODULE_DOC_OUTPUT})
cepgen_generate_python_descriptions(NAME EventModifiers CATEGORY evtmod)
cepgen_generate_python_descriptions(NAME OutputModules CATEGORY evtout)
cepgen_generate_python_descriptions(NAME Integrators CATEGORY integr)
cepgen_generate_python_descriptions(NAME Processes CATEGORY proc)
