import Config.Core as cepgen
import Config.ktProcess_cfi as kt
from Config.PDG_cfi import PDG
from Config.generator_cfi import generator

process = kt.process.clone('pptoff',
    processParameters = cepgen.Parameters(
        pair = PDG.muon,
    ),
    inKinematics = cepgen.Parameters(
        pz = (2562.2, 2562.2),
        partonFluxes = (kt.HeavyIonFlux.PhotonElastic, kt.HeavyIonFlux.PhotonElastic),
        heavyIon1 = (208, 82),
        heavyIon2 = (208, 82),
    ),
    outKinematics = kt.process.outKinematics.clone(
        pt = (5.,),
        energy = (0.,),
        rapidity = (-6., 6.),
        eta = (-10., 10.),
        mx = (1.07, 1000.),
        #--- extra cuts on the p1t(l) and p2t(l) plane
        #ptdiff = (0., 2.5),
        #--- distance in rapidity between l^+ and l^-
        #dely = (4., 5.),
    ),
)

#--- events generation
generator.numEvents = 100000
