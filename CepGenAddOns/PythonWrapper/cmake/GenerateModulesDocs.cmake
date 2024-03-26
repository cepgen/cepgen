macro(cepgen_generate_python_descriptions)
    set(options)
    set(one_val NAME)
    set(multi_vals CATEGORY)
    cmake_parse_arguments(ARG "${options}" "${one_val}" "${multi_vals}" ${ARGN})
    exec_program(${MODULE_DOC_EXE}
        ARGS --documentation-generator "\"text<modulesOnly\""
             --categories ${ARG_CATEGORY}
             --quiet
        OUTPUT_VARIABLE modules)
    foreach(_mod ${modules})
        exec_program(${MODULE_DOC_EXE}
            ARGS --documentation-generator python
                 --categories ${ARG_CATEGORY}
                 --modules ${_mod}
                 --output "${ARGS_NAME}/${_mod}_cfi.py"
                 --quiet
            OUTPUT_VARIABLE dummy)
        message(STATUS "      ... generated Python configuration file for ${ARG_NAME}/${_mod}")
    endforeach()
    message(STATUS "Python... added modules description for ${ARG_CATEGORY} in ${ARG_NAME}")
endmacro()

cepgen_generate_python_descriptions(NAME OutputModules CATEGORY evtout)
cepgen_generate_python_descriptions(NAME Processes CATEGORY proc)

