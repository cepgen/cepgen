import Config.Core as cepgen
from Config.PDG_cfi import PDG

process = cepgen.Module('pptoff',
    processParameters = cepgen.Parameters(
        ktFactorised = True,
        mode = cepgen.ProcessMode.InelasticInelastic,
        pair = PDG.top,
    ),
    inKinematics = cepgen.Parameters(
        cmEnergy = 13.e3,
        structureFunctions = cepgen.StructureFunctions.LUXlike,
    ),
)

