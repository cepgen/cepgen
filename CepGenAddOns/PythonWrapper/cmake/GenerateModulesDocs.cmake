macro(cepgen_generate_python_descriptions)
    set(options)
    set(one_val NAME)
    set(multi_vals CATEGORIES)
    cmake_parse_arguments(ARG "${options}" "${one_val}" "${multi_vals}" ${ARGN})
    exec_program(${MODULE_DOC_EXE}
        ARGS --documentation-generator "\"text<modulesOnly,camelCaseModulesNames\""
             --categories "\"${ARG_CATEGORIES}\""
             --quiet
        OUTPUT_VARIABLE modules)
    message(STATUS "Python... generating modules descriptions for ${ARG_NAME}")
    set(MODULE_DIR ${MODULE_DOC_OUTPUT}/${ARG_NAME})
    file(MAKE_DIRECTORY ${MODULE_DIR})
    file(WRITE ${MODULE_DIR}/__init__.py
      "##\n"
      "# \\defgroup ${ARG_NAME} ${ARG_NAME} modules steering\n"
      "# \\ingroup python ${ARG_NAME}\n\n"
      "# Automatically generated collection of ${ARG_NAME} Python modules definitions"
      " for the steering of '${ARG_CATEGORIES}' category(ies).")
    foreach(_mod ${modules})
        exec_program(${MODULE_DOC_EXE}
            ARGS --documentation-generator python
                 --categories "\"${ARG_CATEGORIES}\""
                 --modules ${_mod}
                 --output "${MODULE_DIR}/${_mod}_cfi.py"
                 --quiet
            OUTPUT_VARIABLE dummy)
        message(DEBUG "      ... generated Python configuration file for ${ARG_NAME}/${_mod}")
    endforeach()
endmacro()

file(MAKE_DIRECTORY ${MODULE_DOC_OUTPUT})
file(WRITE ${MODULE_DOC_OUTPUT}/__init__.py "## \\ingroup python")
cepgen_generate_python_descriptions(NAME Couplings CATEGORIES "alphaem;alphas")
cepgen_generate_python_descriptions(NAME EventModifiers CATEGORIES "evtmod")
cepgen_generate_python_descriptions(NAME FormFactors CATEGORIES "formfac")
cepgen_generate_python_descriptions(NAME Integrators CATEGORIES "integr")
cepgen_generate_python_descriptions(NAME OutputModules CATEGORIES "evtout")
cepgen_generate_python_descriptions(NAME PartonFluxes CATEGORIES "collflux;ktflux")
cepgen_generate_python_descriptions(NAME Processes CATEGORIES "proc")
cepgen_generate_python_descriptions(NAME RandomGenerators CATEGORIES "rndgen")
