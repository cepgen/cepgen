import Config.Core as cepgen
from Config.PDG_cfi import PDG
from Config.timer_cfi import timer # enable timing framework
from Config.generator_cfi import generator as _gen
from Config.ktFluxes_cff import ElectronFlux


process = cepgen.Module('pptoff',
    ktFactorised = True,
    processParameters = cepgen.Parameters(
        pair = PDG.muon,
    ),
    inKinematics = cepgen.Parameters(
        pdgIds = (11, -11),
        pz = (45., 45.),
        partonFluxes = (ElectronFlux.PhotonElasticBudnev, ElectronFlux.PhotonElasticBudnev),
    ),
    outKinematics = cepgen.Parameters(
        pt = (3.,),
    )
)

generator = _gen.clone(  # events generation parameters
    numEvents = 100000,
    printEvery = 10000,
)

text = cepgen.Module('text', # histogramming/ASCII output capability
    histVariables={
        'm(4)': cepgen.Parameters(xrange=(0., 50.), nbins=20),
        'pt(7):pt(8)': cepgen.Parameters(xrange=(0., 25.), yrange=(0., 25.), log=True)
    }
)
dump = cepgen.Module('dump', printEvery = generator.printEvery)
output = cepgen.Sequence(text, dump)
