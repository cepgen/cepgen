import Config.Core as cepgen
from Config.integrators_cff import vegas as integrator
#from Config.pythia8_cff import pythia8 as hadroniser
from Config.ktProcess_cfi import ktProcess

from Config.logger_cfi import logger
logger.enabledModules += ('PPtoFF.prepare',)

process = ktProcess.clone('pptoff',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
    ),
    inKinematics = cepgen.Parameters(
        pz = (6500., 6500.),
        #structureFunctions = cepgen.StructureFunctions.SuriYennie,
        structureFunctions = cepgen.StructureFunctions.LUXlike,
        #structureFunctions = cepgen.StructureFunctions.FioreBrasse,
    ),
    outKinematics = ktProcess.outKinematics.clone(
        pair = 6,
        #eta = (-2.5, 2.5),
        mx = (1.07, 1000.),
        #--- extra cuts on the p1t(t) and p2t(t) plane
        ptdiff = (0., 2000.),
        #--- distance in rapidity between l^+ and l^-
        #dely = (4., 5.),
    ),
)

#--- events generation
from Config.generator_cff import generator
generator.numEvents = 25000
#generator.treat = True

