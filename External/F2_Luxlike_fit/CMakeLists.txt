file(GLOB f2_sources F2_luxlike.f)
file(GLOB as_sources alphaS.f)
file(GLOB mstw_sources mstw*.f)

set(sources ${f2_sources} ${as_sources} ${mstw_sources})
set(ENV{has_schaefer} ${sources})

message(STATUS "LUXlike structure functions found!")

add_library(SchaeferF2 SHARED ${sources})
install(TARGETS SchaeferF2 DESTINATION lib)
