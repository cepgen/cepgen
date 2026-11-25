import Config.Core as cepgen

process = cepgen.Module('pptoww',
    processParameters = cepgen.Parameters(
        ktFactorised = True,
        method = 1,  # legacy steering (now an alias of method 2): must not crash
        mode = cepgen.ProcessMode.ElasticElastic,
    ),
    inKinematics = cepgen.Parameters(
        cmEnergy = 13.e3,
    ),
    #outKinematics = cepgen.Parameters(
    #)
)

