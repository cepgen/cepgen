import Config.Core as cepgen
from Config.Integration.vegas_cff import integrator
from Config.logger_cfi import logger

logger.enabledModules += ('Generator.*',)

import Config.ktProcess_cfi as kt
process = kt.process.clone('pptoww',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.InelasticElastic,
        polarisationStates = 0, # full
        eftExtension = cepgen.Parameters(
            a0w = 1.e-5,
            acw = 0.,
            formFactor = cepgen.Parameters(
                order = -2,
                scale = 500., # GeV
            ),
        ),
    ),
    inKinematics = cepgen.Parameters(
        cmEnergy = 13.e3,
        structureFunctions = cepgen.StructureFunctions.LUXlike,
    ),
    outKinematics = kt.process.outKinematics.clone(
        mx = (1.07, 1000.),
        qt = (0., 1000.),
        #--- extra cuts on the pt(W+) and pt(W-) plane
        ptdiff = (0., 2000.),
        #--- extra cuts on the W+W- system
        invmass = (0.,),
        ptsum = (0.,),
        #--- cuts on single particles' level
        cuts = {
            # cuts on the single W level
            24: cepgen.Parameters(pt = (0.,)), # no pt cut on Ws
            # cuts on the W decay products
            # (mimicking LHC-like experimental cuts)
            #11: cepgen.Parameters(pt = (20.,), eta = (-2.5, 2.5)),
            #13: cepgen.Parameters(pt = (20.,), eta = (-2.5, 2.5))
        },
        #xi = (0.02, 0.15),
    )
)

#--- generation parameters
from Config.generator_cff import generator
generator = generator.clone(
    numEvents = 10000,
    printEvery = 1000,
)
