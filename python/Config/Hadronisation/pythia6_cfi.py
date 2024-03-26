##
# \file
# \defgroup pythia6 Pythia 6 hadronisation module
# \ingroup python pythia6
#
# Base Pythia 6 configuration parameters to define the hadronisation module

from EventModifiers.pythia6_cfi import pythia6 as _pythia6
from Config.Hadronisation.pythia6Tunes_cff import *

pythia6 = _pythia6.clone(
    seed = 1000,
    preConfiguration = (
        'MSTU(21)=1',
    ),
    tuningSettings = pythia6noUESettings,
    #tuningSettings = pythia6CUEP6S1Settings,
    #tuningSettings = pythia6UEZ2starSettings,
    processConfiguration = (
        'tuningSettings',
    ),
)
