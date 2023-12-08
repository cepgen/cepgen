import Config.Core as cepgen
from Config.PDG_cfi import PDG

integrator.numFunctionCalls = 500000
#integrator.verbose = 0

process = cepgen.Module('lpair',
    processParameters = cepgen.Parameters(
        pair = PDG.muon,
    ),
    inKinematics = cepgen.Parameters(
        pz = (15., 15.),
        pdgIds = (11, -11),
    ),
)
