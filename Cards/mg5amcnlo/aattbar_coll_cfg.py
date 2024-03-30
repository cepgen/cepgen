import Config.Core as cepgen
import Config.collinearProcess_cfi as coll
from Config.generator_cfi import generator

#--- process definition
process = coll.process.clone('mg5_aMC',
    processParameters = cepgen.Parameters(
        process = 'a a > t t~',
        mode = cepgen.ProcessMode.ElasticElastic,
    ),
    inKinematics = cepgen.Parameters(
        pz = (6500., 6500.),
        structureFunctions = cepgen.StructureFunctions.luxLike,
    ),
    outKinematics = coll.process.outKinematics.clone(
        q2 = (0., 10.),
        #eta = (-2.5, 2.5),
        mx = (1.07, 1000.),
        pt = (0.,),
    ),
)

#--- events generation
generator.numEvents = 10000

text = cepgen.Module('text',  # histogramming/ASCII output capability
    histVariables={
        'm(4)': cepgen.Parameters(xrange=(0., 250.), nbins=10, log=True),
        'pt(4)': cepgen.Parameters(xrange=(0., 25.), nbins=10, log=True),
        'pt(7):pt(8)': cepgen.Parameters(xrange=(0., 250.), yrange=(0., 250.), log=True)
    }
)
dump = cepgen.Module('dump', printEvery = generator.printEvery)
output = cepgen.Sequence(text, dump)
