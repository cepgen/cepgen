import Config.Core as cepgen
import Config.ktProcess_cfi as ktfactor
from Config.integrators_cff import vegas as integrator
from Config.pdg_cff import PDG

process = ktfactor.process.clone('patoff',
    processParameters = cepgen.Parameters(
        pair = PDG.muon,
    ),
    inKinematics = cepgen.Parameters(
        pz = (6500., 2562.2),
        #structureFunctions = cepgen.StructureFunctions.SuriYennie,
        #structureFunctions = cepgen.StructureFunctions.FioreBrasse,
        #structureFunctions = cepgen.StructureFunctions.ALLM97,
        structureFunctions = cepgen.StructureFunctions.LUXlike,
        ktFluxes = (ktfactor.ProtonFlux.PhotonInelasticBudnev, ktfactor.HeavyIonFlux.PhotonElastic),
        heavyIonB = (208, 82),
    ),
    outKinematics = ktfactor.process.outKinematics.clone(
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

