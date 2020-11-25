import Config.Core as cepgen

#--- process definition
process = cepgen.Module('mg5_aMC',
    processParameters = cepgen.Parameters(
        process = 'a a > mu+ mu-',
        mode = cepgen.ProcessMode.ElasticElastic,
    ),
    inKinematics = cepgen.Parameters(
        pz = (6500., 6500.),
        structureFunctions = cepgen.StructureFunctions.LUXlike,
    ),
    outKinematics = cepgen.Parameters(
        #eta = (-2.5, 2.5),
        mx = (1.07, 2000.),
        ptsingle = (20.,),
    ),
)

#--- events generation
from Config.generator_cff import generator
generator.numEvents = 5000

