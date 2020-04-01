from Config.containers_cfi import Module
from pythia6Tunes_cfi import *

pythia6 = Module('pythia6',
    seed = 1000,
    maxTrials = 1,
    preConfiguration = (
        'MSTU(21)=1',
    ),
    processConfiguration = (
        #'pythia6CUEP6S1Settings',
        'pythia6UEZ2starSettings',
    ),

)
