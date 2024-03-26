import Config.Core as cepgen
from Config.PDG_cfi import PDG
from Config.timer_cfi import timer # enable timing framework
from Config.generator_cfi import generator as _gen
from OutputModules.dump_cfi import dump as _dump_output # periodic event printout
from OutputModules.text_cfi import text as _text_output # ASCII histograms
#from OutputModules.rootTree_cfi import rootTree # dump everything into a flat tree


process = cepgen.Module('lpair',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
        pair = PDG.muon,
    ),
    inKinematics = cepgen.Parameters(
        heavyIon1 = (208, 82),
        heavyIon2 = (208, 82),
        pz = (2562.2, 2562.2),
    ),
    outKinematics = cepgen.Parameters(
        pt = (5.,),
        eta = (-10., 10.)
    )
)

# events generation parameters
generator = _gen.clone(
    numEvents = 100000,
    printEvery = 10000,
)

text = _text_output.clone(
    histVariables={
        'm(4)': cepgen.Parameters(xrange=(5., 25.), nbins=20),
        'pt(7):pt(8)': cepgen.Parameters(xrange=(0., 20.), yrange=(0., 20.), log=True)
    }
)
dump = _dump_output.clone(printEvery = generator.printEvery)
output = cepgen.Sequence(text, dump)
