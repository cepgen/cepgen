import Config.Core as cepgen
from Config.integrators_cff import vegas as integrator
from Config.ktProcess_cfi import ktProcess, KTFlux

process = ktProcess.clone('patoll',
    inKinematics = cepgen.Parameters(
        pz = (6500., 2562.2),
        structureFunctions = cepgen.StructureFunctions.SuriYennie,
        #structureFunctions = cepgen.StructureFunctions.FioreBrasse,
        #structureFunctions = cepgen.StructureFunctions.ALLM91,
        ktFluxes = (KTFlux.PhotonInelasticBudnev, KTFlux.PhotonElasticHI),
        heavyIonB = (208, 82),
    ),
    outKinematics = ktProcess.outKinematics.clone(
        pair = 13,
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

#--- events generation
from Config.generator_cff import generator
generator.numEvents = 100000

