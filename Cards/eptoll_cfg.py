import Config.Core as cepgen
import Config.ktProcess_cfi as kt
from Config.PDG_cfi import PDG
from Config.generator_cfi import generator
from FormFactors.standardDipole_cfi import standardDipole
from FormFactors.pointLikeFermion_cfi import pointLikeFermion

process = kt.process.clone('pptoff',
    processParameters = cepgen.Parameters(
        pair = PDG.muon,
    ),
    inKinematics = cepgen.Parameters(
        pdgIds = (PDG.proton, PDG.electron),
        formFactors = [standardDipole, pointLikeFermion],
        pz = (7000., 50.),
        structureFunctions = cepgen.StructureFunctions.suriYennie,
    ),
    outKinematics = kt.process.outKinematics.clone(
        pt = (2.5,),
        energy = (0.,),
        eta = (-2.5, 2.5),
        mx = (1.07, 1000.),
    ),
)

generator.numEvents = 10000
text = cepgen.Module('text', # histogramming/ASCII output capability
    histVariables={
        'm(4)': cepgen.Parameters(xrange=(0., 25.), nbins=20, log=False),
        'pt(7):pt(8)': cepgen.Parameters(xrange=(0., 50.), yrange=(0., 50.), log=True)
    }
)
output = cepgen.Sequence(text)
