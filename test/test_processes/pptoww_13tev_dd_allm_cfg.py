import Config.Core as cepgen
from Config.Integration.vegas_cff import integrator

process = cepgen.Module('pptoww',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.InelasticInelastic,
    ),
    inKinematics = cepgen.Parameters(
        cmEnergy = 13.e3,
        structureFunctions = cepgen.StructureFunctions.ALLM97,
    ),
    outKinematics = cepgen.Parameters(
        mx = (1.07, 1000.),
    )
)

