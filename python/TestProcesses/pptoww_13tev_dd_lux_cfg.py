import Config.Core as cepgen

process = cepgen.Module('pptoww',
    processParameters = cepgen.Parameters(
        ktFactorised = True,
        mode = cepgen.ProcessMode.InelasticInelastic,
    ),
    inKinematics = cepgen.Parameters(
        cmEnergy = 13.e3,
        structureFunctions = cepgen.StructureFunctions.LUXlike,
    ),
    outKinematics = cepgen.Parameters(
        mx = (1.07, 1000.),
    )
)

