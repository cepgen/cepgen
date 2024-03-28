cepgen_test_category(NAME "Validation" SOURCES *.cc)

add_test(NAME cross_section_vegas COMMAND processes_cross_sections -i Vegas -q --cfg test_processes.cfg)
#cepgen_test_declare_unstable(processes_cross_sections cross_section_vegas)  # FIXME remove once agreement is reached
