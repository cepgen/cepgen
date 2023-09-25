import Config.Core as cepgen
import Config.collinearProcess_cfi as coll
from Config.generator_cff import generator
from Config.Logger_cfi import logger
#logger.enabledModules += ('MadGraphProcess.eval',)

#--- process definition
process = coll.process.clone('mg5_aMC',
    processParameters = cepgen.Parameters(
        process = 'a a > mu+ mu-',
        # alternatively, if shared object is already generated
        #lib = 'libCepGenMadGraphProcess.so',
        # alternatively, if standalone_cpp directory is already generated
        #standaloneCppPath = '/tmp/cepgen_mg5_aMC',
        mode = cepgen.ProcessMode.ElasticElastic,
    ),
    inKinematics = cepgen.Parameters(
        pz = (6500., 6500.),
        structureFunctions = cepgen.StructureFunctions.LUXlike,
    ),
    outKinematics = coll.process.outKinematics.clone(
        #eta = (-2.5, 2.5),
        qt = (0., 10.),
        mx = (1.07, 1000.),
        pt = (0.,),
    ),
)

#--- events generation
generator.numEvents = 10000

text = cepgen.Module('text',  # histogramming/ASCII output capability
    #variables = ['nev', 'm(4)', 'tgen'],
    histVariables={
        'm(4)': cepgen.Parameters(xrange=(0., 250.), nbins=10, log=True),
        'pt(4)': cepgen.Parameters(xrange=(0., 25.), nbins=10, log=True),
        'pt(7):pt(8)': cepgen.Parameters(xrange=(0., 250.), yrange=(0., 250.), log=True)
    }
)
dump = cepgen.Module('dump', printEvery = generator.printEvery)
output = cepgen.Sequence(text, dump)
