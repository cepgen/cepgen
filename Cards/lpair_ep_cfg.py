import Config.Core as cepgen
from Config.PDG_cfi import PDG
from Config.timer_cfi import timer # enable timing framework
from Config.generator_cfi import generator as _gen
from FormFactors.standardDipole_cfi import standardDipole
from FormFactors.pointLikeFermion_cfi import pointLikeFermion
#from OutputModules.rootTree_cfi import rootTree # dump everything into a flat tree


process = cepgen.Module('lpair',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
        pair = PDG.muon,
    ),
    inKinematics = cepgen.Parameters(
        pdgIds = (2212, 11),
        formFactors = [standardDipole, pointLikeFermion],
        pz = (7000., 50.),
    ),
    outKinematics = cepgen.Parameters(
        pt = (2.5,),
    )
)

# events generation parameters
generator = _gen.clone(
    numEvents = 10000,
    printEvery = 2500,
)

text = cepgen.Module('text', # histogramming/ASCII output capability
    histVariables={
        'm(4)': cepgen.Parameters(xrange=(0., 25.), nbins=20, log=False),
        'pt(7):pt(8)': cepgen.Parameters(xrange=(0., 5.), yrange=(0., 5.), log=True)
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
