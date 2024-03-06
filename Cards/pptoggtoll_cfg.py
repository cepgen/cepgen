import Config.Core as cepgen
import Config.ktProcess_cfi as kt
from Config.PDG_cfi import PDG
from Config.generator_cfi import generator

process = kt.process.clone('pptoff',
    processParameters = cepgen.Parameters(
        pair = PDG.muon,
        mode = cepgen.ProcessMode.InelasticInelastic,
    ),
    inKinematics = cepgen.Parameters(
        pdgIds = (PDG.proton, PDG.proton),
        pz = (6500., 6500.),
        partonFluxes = kt.ProtonFlux.GluonKMR
    ),
    outKinematics = kt.process.outKinematics.clone(
        pt = (25.,),
        energy = (0.,),
        eta = (-2.5, 2.5),
        mx = (1.07, 1000.),
        #--- extra cuts on the p1t(l) and p2t(l) plane
        #ptdiff = (0., 2.5),
        #--- distance in rapidity between l^+ and l^-
        #dely = (4., 5.),
    ),
)

generator.numEvents = 100000

root = cepgen.Module('root_tree')
dump = cepgen.Module('dump', printEvery = generator.printEvery)
output = cepgen.Sequence(
    root,
    dump,
)
