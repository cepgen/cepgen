import Config.Core as cepgen
import Config.collinearProcess_cfi as coll
from Config.generator_cfi import generator
from FormFactors.standardDipole_cfi import standardDipole
from FormFactors.pointLikeFermion_cfi import pointLikeFermion

#--- process definition
process = coll.process.clone('mg5_aMC:eh',
    processParameters = cepgen.Parameters(
        process = 'a e- > ta+ ta- e- / h z e-',
        model = 'sm-full',
        mode = cepgen.ProcessMode.ElasticElastic,
        #kinematicsGenerator = cepgen.Module('kt_single:2to4'),
        kinematicsGenerator = cepgen.Module('kt:2to4'),
    ),
    inKinematics = cepgen.Parameters(
        pdgIds = (2212, 11),
        formFactors = [standardDipole, pointLikeFermion],
        pz = (7000., 50.),
        structureFunctions = cepgen.StructureFunctions.luxLike,
    ),
    outKinematics = coll.process.outKinematics.clone(
        q2 = (0., 10.),
        #eta = (-2.5, 2.5),
        mx = (1.07, 1000.),
        pt = (2.5,),
    ),
)

generator.numEvents = 10000
dump = cepgen.Module('dump', printEvery = generator.printEvery)
output = cepgen.Sequence(dump)
