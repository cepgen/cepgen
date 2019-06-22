from Config.containers_cfi import Module, Parameters

pythia6 = Module('pythia6',
    seed = 1000,
    maxTrials = 1,
    preConfiguration = (
        'MSTU(21)=1',
    ),
    processConfiguration = (
        #'commonSettings', 'tuningSettings',
    ),

)
