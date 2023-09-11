import Config.Core as cepgen
from Config.PDG_cfi import PDG

process = cepgen.Module('pptoff',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
        pair = PDG.top,
    ),
    inKinematics = cepgen.Parameters(
        cmEnergy = 13.e3,
    ),
)

