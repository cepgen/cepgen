import Config.Core as cepgen
import Config.ktProcess_cfi as kt
from Config.PDG_cfi import PDG
from Config.generator_cfi import generator

#--- process definition
process = kt.process.clone('mg5_aMC',
    processParameters = cepgen.Parameters(
        process = 'a a > w+ w-',
        model = 'QAll_5_Aug21v2',
        # alternatively, if shared object is already generated
        #lib = 'libCepGenMadGraphProcess.so',
        # alternatively, if standalone_cpp directory is already generated
        #standaloneCppPath = '/tmp/cepgen_mg5_aMC',
        mode = cepgen.ProcessMode.ElasticElastic,
        modelParameters = cepgen.Parameters(
            anoinputs = {4: 1.e-4}
        )
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
