##
# \file
# \ingroup python pythia6
#
# Base Pythia 6 configuration parameters to define the hadronisation module


from Config.containers_cff import Module
from Config.Hadronisation.pythia6Tunes_cff import *

pythia6 = Module('pythia6',
    seed = 1000,
    maxTrials = 1,
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
