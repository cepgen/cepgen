import Config.Core as cepgen
from Config.Integration.vegas_cff import integrator

process = cepgen.Module('lpair',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.InelasticElastic,
        pair = 13,
    ),
    inKinematics = cepgen.Parameters(
        cmEnergy = 13.e3,
        structureFunctions = cepgen.StructureFunctions.SuriYennie,
    ),
    outKinematics = cepgen.Parameters(
        pt = (15.,),
        eta = (-2.5, 2.5),
        mx = (1.07, 1000.),
    )
)

