import Config.Core as cepgen
import Config.ktProcess_cfi as kt
from Config.PDG_cfi import PDG, registerParticle
from Config.generator_cfi import generator
from OutputModules.dump_cfi import dump as _dump_output # periodic event printout
from OutputModules.text_cfi import text as _text_output # ASCII histograms

#--- two-photon production of a chargino pair, gamma gamma -> chi_1^+ chi_1^-
# A chargino is a charged Dirac fermion: its photon coupling only depends on its charge and mass, i.e. the process is
# the same as gamma gamma -> f fbar, for a (colourless) fermion of charge +-1 and mass m(chi_1^+).
registerParticle(1000024, 'chargino_1', mass=150., charge=1., colour=0, fermion=True)

process = kt.process.clone('pptoff',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
        pair = PDG.chargino_1,
        method = 0, # on-shell matrix element
    ),
    inKinematics = cepgen.Parameters(
        pdgIds = (PDG.proton, PDG.proton),
        pz = (6500., 6500.),
    ),
    outKinematics = kt.process.outKinematics.clone(
        pt = (0.,),
        eta = (-2.5, 2.5),
    ),
)

#--- events generation
generator.numEvents = 10000
text = _text_output.clone(
    histVariables={
        'm(4)': cepgen.Parameters(xrange=(300., 1000.), nbins=14),
    }
)
dump = _dump_output.clone(printEvery = generator.printEvery)
output = cepgen.Sequence(text, dump)
