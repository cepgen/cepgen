import Config.Core as cepgen
from Config.Integration.vegas_cff import integrator
from Config.logger_cfi import logger

'''from Config.Hadronisation.pythia8_cff import pythia8
hadroniser = pythia8.clone('pythia8',
    preConfiguration = pythia8.preConfiguration+(
        #'PartonLevel:MPI = on',
        #'PartonLevel:ISR = on',
        #'PartonLevel:FSR = on',
        'ProcessLevel:resonanceDecays = off', # disable the W decays
    ),
    pythiaConfiguration = (
        # process-specific
        '13:onMode = off', # disable muon decays
        '24:onMode = off', # disable all W decays, but...
        #'24:onIfAny = 11 13', # enable e-nue + mu-numu final states
        '24:onPosIfAny = 11', # enable W- -> e- + nu_e decay
        '24:onNegIfAny = 13', # enable W+ -> mu+ + nu_mu decay
    ),
    processConfiguration = pythia8.processConfiguration+('pythiaConfiguration',),
)'''
logger.enabledModules += ('Hadroniser.configure', 'Generator.*',)

import Config.ktProcess_cfi as kt
process = kt.process.clone('pptoww',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.InelasticElastic,
        polarisationStates = 0, # full
    ),
    inKinematics = cepgen.Parameters(
        cmEnergy = 13.e3,
        #structureFunctions = cepgen.StructureFunctions.SzczurekUleshchenko,
        #structureFunctions = cepgen.StructureFunctions.ALLM97,
        structureFunctions = cepgen.StructureFunctions.LUXlike,
    ),
    outKinematics = kt.process.outKinematics.clone(
        mx = (1.07, 1000.),
        qt = (0., 1000.),
        #--- extra cuts on the pt(W+) and pt(W-) plane
        ptdiff = (0., 2000.),
        #--- extra cuts on the W+W- system
        invmass = (0.,),
        ptsum = (0.,),
        #--- cuts on single particles' level
        cuts = {
            # cuts on the single W level
            24: cepgen.Parameters(pt = (0.,)), # no pt cut on Ws
            # cuts on the W decay products
            # (mimicking LHC-like experimental cuts)
            11: cepgen.Parameters(pt = (20.,), eta = (-2.5, 2.5)),
            13: cepgen.Parameters(pt = (20.,), eta = (-2.5, 2.5))
        },
        #xi = (0.02, 0.15),
    )
)

#--- generation parameters
from Config.generator_cff import generator
generator = generator.clone(
    numEvents = 10000,
    printEvery = 1000,
)
