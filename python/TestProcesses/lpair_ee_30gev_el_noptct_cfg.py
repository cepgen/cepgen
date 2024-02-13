import Config.Core as cepgen
from Config.PDG_cfi import PDG

process = cepgen.Module('lpair',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
        pair = PDG.muon,
    ),
    inKinematics = cepgen.Parameters(
        pz = (15., 15.),
        pdgIds = (11, -11),
    ),
)
