file(GLOB sources F2_luxlike.f alphaS.f mstw*.f)

set(ENV{has_schaefer} ${sources})
message(STATUS "LUXlike structure functions found!")

