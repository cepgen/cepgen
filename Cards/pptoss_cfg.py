import Config.Core as cepgen
import Config.ktProcess_cfi as kt
from Config.PDG_cfi import PDG, registerParticle
from Config.generator_cfi import generator
from OutputModules.dump_cfi import dump as _dump_output # periodic event printout
from OutputModules.text_cfi import text as _text_output # ASCII histograms

#--- two-photon production of a pair of charged scalars
registerParticle(2000013, 'smuon_R', mass=100., charge=1., colour=0, fermion=False)

process = kt.process.clone('pptoss',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
        pair = PDG.smuon_R,
    ),
    inKinematics = cepgen.Parameters(
        pdgIds = (PDG.proton, PDG.proton),
        sqrtS = 14.e3,
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
        'm(4)': cepgen.Parameters(xrange=(200., 1000.), nbins=16),
        'pt(7):pt(8)': cepgen.Parameters(xrange=(0., 250.), yrange=(0., 250.), log=True)
    }
)
dump = _dump_output.clone(printEvery = generator.printEvery)
output = cepgen.Sequence(text, dump)
