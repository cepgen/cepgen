import Config.Core as cepgen
from Config.integrators_cff import vegas as integrator
from Config.ktProcess_cfi import ktProcess
from Config.logger_cfi import logger
from Config.pythia8_cff import pythia8 as hadroniser

hadroniser.pythiaProcessConfiguration += (
    # process-specific
    '13:onMode = off', # disable muon decays
    '24:onMode = off', # disable all W decays, but...
    #'24:onIfAny = 11 13', # enable e-nue + mu-numu final states
    '24:onPosIfAny = 11',
    '24:onNegIfAny = 13', # enable e-nue + mu-numu final states
)
hadroniser.pythiaPreConfiguration += (
    #'PartonLevel:MPI = on',
    #'PartonLevel:ISR = on',
    'PartonLevel:FSR = on',
    'ProcessLevel:resonanceDecays = off',
#    'BeamRemnants:unresolvedHadron = 3',
#    'Photon:ProcessType = 4',
    #'PartonLevel:Remnants = off',
#    'TimeShower:MEcorrections = off',
#    'TimeShower:globalRecoil = on',
)

#logger.level = cepgen.Logging.DebugInsideLoop

process = ktProcess.clone('pptoww',
    mode = cepgen.ProcessMode.InelasticInelastic,
    inKinematics = cepgen.Parameters(
        cmEnergy = 13.e3,
        #structureFunctions = cepgen.StructureFunctions.SuriYennie,
        #structureFunctions = cepgen.StructureFunctions.FioreBrasse,
        #structureFunctions = cepgen.StructureFunctions.SzczurekUleshchenko,
        #structureFunctions = cepgen.StructureFunctions.ALLM91,
        #structureFunctions = cepgen.StructureFunctions.ALLM97,
        structureFunctions = cepgen.StructureFunctions.LUXlike,
    ),
    outKinematics = ktProcess.outKinematics.clone(
        mx = (1.07, 1000.),
        qt = (0., 1000.),
        #--- extra cuts on the pt(W+) and pt(W-) plane
        ptdiff = (0., 2000.),
        #--- distance in rapidity between W+ and W-
        #rapiditydiff = (4., 5.),
        cuts = {
            # cuts on the single W level
            24: cepgen.Parameters(
                pt = (0.,),
            ),
            # cuts on the W decay products
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
from Config.generator_cff import generator
generator.numEvents = 1000
generator.printEvery = 100
