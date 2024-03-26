import Config.Core as cepgen
import Config.ktProcess_cfi as kt
from Config.PDG_cfi import PDG
from Config.generator_cfi import generator
from Integrators.vegas_cfi import vegas as integrator


process = kt.process.clone('pptoff',
    processParameters = cepgen.Parameters(
        pair = PDG.muon,
    ),
    inKinematics = cepgen.Parameters(
        beam1id = 2212,
        heavyIon2 = (208, 82),
        pz = (6500., 2562.2),
        #structureFunctions = cepgen.StructureFunctions.SuriYennie,
        #structureFunctions = cepgen.StructureFunctions.FioreBrasse,
        #structureFunctions = cepgen.StructureFunctions.ALLM97,
        structureFunctions = cepgen.StructureFunctions.LUXlike,
        partonFluxes = (kt.ProtonFlux.PhotonInelasticBudnev, kt.HeavyIonFlux.PhotonElastic),
    ),
    outKinematics = kt.process.outKinematics.clone(
        pt = (4.,),
        energy = (0.,),
        rapidity = (-6., 7.),
        #eta = (-2.5, 2.5),
        mx = (1.07, 1000.),
        #--- extra cuts on the p1t(l) and p2t(l) plane
        #ptdiff = (0., 2.5),
        #--- distance in rapidity between l^+ and l^-
        #dely = (4., 5.),
    ),
)

generator.numEvents = 100000
