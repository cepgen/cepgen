import Config.Core as cepgen
from Config.Integration.vegas_cff import integrator

process = cepgen.Module('pptoww',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
    ),
    inKinematics = cepgen.Parameters(
        cmEnergy = 13.e3,
    ),
    #outKinematics = cepgen.Parameters(
    #)
)

