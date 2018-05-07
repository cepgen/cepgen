import Config.Core as cepgen
from Config.integrators_cff import vegas as integrator
from Config.ktProcess_cfi import ktProcess
from Config.logger_cfi import logger
from Config.pythia8_cff import pythia8 as hadroniser
from Config.generator_cff import generator

hadroniser.pythiaProcessConfiguration += (
    # process-specific
    '13:onMode = off', # disable muon decays
    '24:onMode = off', # disable all W decays, but...
    #'24:onIfAny = 11 13', # enable e-nue + mu-numu final states
    '24:onPosIfAny = 11', # enable W- -> e- + nu_e decay
    '24:onNegIfAny = 13', # enable W+ -> mu+ + nu_mu decay
)
hadroniser.pythiaPreConfiguration += (
    #'PartonLevel:MPI = on',
    #'PartonLevel:ISR = on',
    #'PartonLevel:FSR = on',
    'ProcessLevel:resonanceDecays = off', # disable the W decays
)

#logger.level = cepgen.Logging.DebugInsideLoop

process = ktProcess.clone('pptoww',
    mode = cepgen.ProcessMode.InelasticInelastic,
    inKinematics = cepgen.Parameters(
        cmEnergy = 13.e3,
        #structureFunctions = cepgen.StructureFunctions.SzczurekUleshchenko,
        #structureFunctions = cepgen.StructureFunctions.ALLM97,
        structureFunctions = cepgen.StructureFunctions.LUXlike,
    ),
    outKinematics = ktProcess.outKinematics.clone(
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
            24: cepgen.Parameters(
                pt = (0.,), # no pt cut on Ws
            ),
            # cuts on the W decay products
            # (mimicking LHC-like experimental cuts)
            11: cepgen.Parameters(
                pt = (20.,),
                eta = (-2.5, 2.5),
            ),
            13: cepgen.Parameters(
                pt = (20.,),
                eta = (-2.5, 2.5),
            )
        },
    )
)

#--- import the default generation parameters
generator = generator.clone(
    numEvents = 1000,
    printEvery = 100,
    numThreads = 4,
    treat = True, # smoothing of the integrand
)
