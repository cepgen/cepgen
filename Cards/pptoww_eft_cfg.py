import Config.Core as cepgen
import Config.ktProcess_cfi as kt
from Config.PDG_cfi import PDG
from Config.Integration.vegas_cfi import vegas as integrator
from Config.logger_cfi import logger
from Config.generator_cfi import generator as _gen

#logger.enabledModules += ('Generator.*',)


class EFTModel:
    SM = 0
    W = 1
    Wbar = 2
    phiW = 3
    phiWbar = 4
    phiB = 5
    phiBbar = 6
    WB = 7
    WbarB = 8


process = kt.process.clone('pptoww',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
        polarisationStates = 0, # full
        model = EFTModel.W,
        eftParameters = cepgen.Parameters(
            s1 = 1.e-2
        ),
    ),
    inKinematics = cepgen.Parameters(
        pdgIds = (PDG.proton, PDG.proton),
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

generator = _gen.clone(  # generation parameters
    numEvents = 10000,
    printEvery = 1000,
)
