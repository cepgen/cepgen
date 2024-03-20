import Config.Core as cepgen
from Config.PDG_cfi import PDG
from Config.timer_cfi import timer # enable timing framework
from Config.generator_cfi import generator as _gen
#from Config.OutputModule.rootTree_cfi import rootTree # dump everything into a flat tree


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

text = cepgen.Module('text', # histogramming/ASCII output capability
    histVariables={
        'm(4)': cepgen.Parameters(xrange=(5., 25.), nbins=20),
        'pt(7):pt(8)': cepgen.Parameters(xrange=(0., 20.), yrange=(0., 20.), log=True)
    }
)
#lhef = cepgen.Module('lhef', filename='test.lhe')
#hepmc = cepgen.Module('hepmc', filename='test.hepmc')
dump = cepgen.Module('dump', printEvery = generator.printEvery)
output = cepgen.Sequence(
    #rootTree,
    text,
    #lhef,
    #hepmc,
    dump,
)
