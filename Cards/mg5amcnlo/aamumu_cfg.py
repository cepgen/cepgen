import Config.Core as cepgen
import Config.ktProcess_cfi as kt
from Config.PDG_cfi import PDG
from Config.generator_cfi import generator
#from Config.logger_cfi import logger
#logger.enabledModules += ('MadGraphProcess.eval',)

#--- process definition
process = kt.process.clone('mg5_aMC',
    processParameters = cepgen.Parameters(
        process = 'a a > mu+ mu-',
        # alternatively, if shared object is already generated
        #lib = 'libCepGenMadGraphProcess.so',
        # alternatively, if standalone_cpp directory is already generated
        #standaloneCppPath = '/tmp/cepgen_mg5_aMC',
        mode = cepgen.ProcessMode.ElasticElastic,
    ),
    inKinematics = cepgen.Parameters(
        pdgIds = (PDG.proton, PDG.proton),
        pz = (6500., 6500.),
        structureFunctions = cepgen.StructureFunctions.luxLike,
    ),
    outKinematics = kt.process.outKinematics.clone(
        qt = (0., 10.),
        #eta = (-2.5, 2.5),
        mx = (1.07, 1000.),
        pt = (25.,),
    ),
)

#--- events generation
generator.numEvents = 10000

dump = cepgen.Module('dump', printEvery = generator.printEvery)
output = cepgen.Sequence(dump)
