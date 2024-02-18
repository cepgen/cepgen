import Config.Core as cepgen

process = cepgen.Module('pptoff',
    processParameters = cepgen.Parameters(
        ktFactorised = True,
        mode = cepgen.ProcessMode.ElasticElastic,
        pair = 13,
        offShellParameters = cepgen.Parameters(mat1 = 2, mat2 = 0)
    ),
    inKinematics = cepgen.Parameters(
        cmEnergy = 13.e3,
    ),
    outKinematics = cepgen.Parameters(
        pt = (3.,),
        eta = (-2.5, 2.5),
    )
)

