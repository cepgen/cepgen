import Config.Core as cepgen
import Config.ktProcess_cfi as kt
from Config.integrators_cff import vegas as integrator
from Config.PDG_cfi import PDG

from Config.logger_cfi import logger
logger.enabledModules += ('GenericKTProcess.registerVariable',)

process = kt.process.clone('patoff',
    processParameters = cepgen.Parameters(
        pair = PDG.charm,
    ),
    inKinematics = cepgen.Parameters(
        pz = (6500., 2562.2),
        structureFunctions = cepgen.StructureFunctions.SuriYennie,
        #structureFunctions = cepgen.StructureFunctions.FioreBrasse,
        ktFluxes = (kt.ProtonFlux.GluonKMR, kt.HeavyIonFlux.PhotonElastic),
        #ktFluxes = (kt.ProtonFlux.PhotonElastic, kt.HeavyIonFlux.PhotonElastic),
        heavyIonB = (208, 82),
        kmrGridPath = 'gluon_mmht2014nlo_Watt.dat',
    ),
    outKinematics = kt.process.outKinematics.clone(
        pt = (0.,),
        energy = (0.,),
        rapidity = (-7., 9.),
        #qt = (0.,1000.),
        #eta = (-2.5, 2.5),
        mx = (1.07, 1000.),
        #--- extra cuts on the p1t(l) and p2t(l) plane
        #ptdiff = (0., 2.5),
        #--- distance in rapidity between l^+ and l^-
        #dely = (4., 5.),
    ),
)

#--- events generation
from Config.generator_cff import generator
generator.numEvents = 10000
generator.numThreads = 1

