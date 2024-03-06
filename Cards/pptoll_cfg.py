import Config.Core as cepgen
import Config.ktProcess_cfi as kt
#from Config.Hadronisation.pythia6_cfi import pythia6 as hadroniser
#from Config.Hadronisation.pythia8_cfi import pythia8 as hadroniser
from Config.PDG_cfi import PDG
from Config.generator_cfi import generator

#--- example of an auxiliary particles definition
#registerParticle(1000001, 'sd_l', mass=100., charge=1., fermion=True) # right now, only fermionic coupling handled

process = kt.process.clone('pptoff',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.InelasticElastic,
        pair = PDG.muon,
        #pair = PDG.up,
        #pair = PDG.sd_l, # whatever was defined above as "new" particle
    ),
    inKinematics = cepgen.Parameters(
        pdgIds = (PDG.proton, PDG.proton),
        #pdgIds = (PDG.proton, PDG.electron),
        pz = (6500., 6500.),
        structureFunctions = cepgen.StructureFunctions.SuriYennie,
        #structureFunctions = cepgen.StructureFunctions.FioreBrasse,
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

#--- events generation
generator.numEvents = 10000
