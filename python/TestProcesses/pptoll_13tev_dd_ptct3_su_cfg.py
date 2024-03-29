import Config.Core as cepgen

process = cepgen.Module('pptoff',
    processParameters = cepgen.Parameters(
        ktFactorised = True,
        mode = cepgen.ProcessMode.InelasticInelastic,
        pair = 13,
        offShellParameters = cepgen.Parameters(mat1 = 2, mat2 = 0)
    ),
    inKinematics = cepgen.Parameters(
        cmEnergy = 13.e3,
        structureFunctions = cepgen.StructureFunctions.SzczurekUleshchenko,
    ),
    outKinematics = cepgen.Parameters(
        pt = (3.,),
        eta = (-2.5, 2.5),
        mx = (1.07, 1000.),
    )
)

