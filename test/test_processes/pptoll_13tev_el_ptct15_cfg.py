import Config.Core as cepgen
from Config.Integration.vegas_cff import integrator

process = cepgen.Module('pptoff',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
        pair = 13,
    ),
    inKinematics = cepgen.Parameters(
        cmEnergy = 13.e3,
    ),
    outKinematics = cepgen.Parameters(
        pt = (15.,),
        eta = (-2.5, 2.5),
    )
)

