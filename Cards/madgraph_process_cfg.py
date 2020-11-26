import Config.Core as cepgen
import Config.ktProcess_cfi as kt

#--- process definition
#process = cepgen.Module('mg5_aMC',
process = kt.process.clone('mg5_aMC',
    processParameters = cepgen.Parameters(
        process = 'a a > mu+ mu-',
        # alternatively, if shared object is already generated
        #lib = 'libCepGenMadGraphProcess.so',
        # alternatively, if standalone_cpp directory is already generated
        standaloneCppPath = '/tmp/cepgen_mg5_aMC',
        mode = cepgen.ProcessMode.ElasticElastic,
    ),
    inKinematics = cepgen.Parameters(
        pz = (6500., 6500.),
        structureFunctions = cepgen.StructureFunctions.LUXlike,
    ),
    #outKinematics = cepgen.Parameters(
    outKinematics = kt.process.outKinematics.clone(
        #eta = (-2.5, 2.5),
        qt = (0., 10.),
        mx = (1.07, 1000.),
        pt = (0.,),
    ),
)

print(process)

#--- events generation
from Config.generator_cff import generator
generator.numEvents = 5000

