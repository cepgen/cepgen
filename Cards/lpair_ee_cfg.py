import Config.Core as cepgen
from Config.PDG_cfi import PDG
from Config.Timer_cfi import timer # enable timing framework
from Config.Cuts_cfi import vermaserenCuts

process = cepgen.Module('lpair',
    processParameters = cepgen.Parameters(
        pair = PDG.muon,
    ),
    inKinematics = cepgen.Parameters(
        pz = (45., 45.),
        pdgIds = (11, -11),
    ),
    outKinematics = vermaserenCuts,
)

# events generation parameters
from Config.generator_cff import generator
generator = generator.clone(
    numEvents = 100000,
    printEvery = 10000,
)

#from Config.OutputModule.ROOTTree_cfi import rootTree # dump everything into a flat tree
text = cepgen.Module('text', # histogramming/ASCII output capability
    histVariables={
        'm(4)': cepgen.Parameters(xrange=(0., 10.), nbins=20, log=True),
        'pt(7):pt(8)': cepgen.Parameters(xrange=(0., 10.), yrange=(0., 10.), log=True)
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
