##
# \file
# \defgroup herwig6 Herwig 6 hadronisation module
# \ingroup python herwig6
#
# Base Herwig 6 configuration parameters to define the hadronisation module

from EventModifiers.herwig6_cfi import herwig6 as _herwig6

herwig6 = _herwig6.clone(
    seed = 1000,
    # any parameter steered here is applied before the first event is hadronised,
    # e.g. 'BTCLM=1.0', 'CLMAX=3.35', 'CLPOW=2.0', 'PSPLT(1)=1.0', 'IOPREM=1'
    preConfiguration = (),
    processConfiguration = (),
)
